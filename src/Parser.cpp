
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

#include "Parser.h"
#include "Filament.h"
#include "Cylinder.h"
#include "Bead.h"

#include "SysParams.h"

namespace medyan {

void SystemParser::initInputHeader() {
    headerParser.addComment("##################################################");
    headerParser.addComment("### Important notes:");
    headerParser.addComment("### 1. Units in MEDYAN are nm, second, pN, and pN*nm");
    headerParser.addComment("##################################################");
    headerParser.addEmptyLine();
}

void SystemParser::initChemParser() {
    using namespace std;

    chemParser.addComment("##################################################");
    chemParser.addComment("### Chemistry parameters");
    chemParser.addComment("##################################################");
    chemParser.addEmptyLine();

    chemParser.addComment("###### Chemistry algorithm ######");
    chemParser.addEmptyLine();

    chemParser.addStringArgsWithAliases(
        "CALGORITHM", { "CALGORITHM:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            if(lineVector.size() != 2) {
                cout <<
                     "There was an error parsing input file at Chemistry algorithm. Exiting."
                     << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                sc.chemParams.chemistryAlgorithm.algorithm = lineVector[1];
            }
        },
        [] (const SimulConfig& sc) {
            return vector<string> { sc.chemParams.chemistryAlgorithm.algorithm };
        }
    );
    chemParser.addEmptyLine();

    chemParser.addComment("# Use either time mode or step mode");
    chemParser.addStringArgsWithAliases(
        "RUNTIME", { "RUNTIME:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            if(lineVector.size() > 2) {
                cout <<
                     "There was an error parsing input file at Chemistry algorithm. Exiting."
                     << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                sc.chemParams.chemistryAlgorithm.runTime = atof(lineVector[1].c_str());
            }
        },
        [] (const SimulConfig& sc) {
            return vector<string> { to_string(sc.chemParams.chemistryAlgorithm.runTime) };
        }
    );
    chemParser.addStringArgsWithAliases(
        "SNAPSHOTTIME", { "SNAPSHOTTIME:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            if(lineVector.size() > 2) {
                cout <<
                     "There was an error parsing input file at Chemistry algorithm. Exiting."
                     << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                sc.chemParams.chemistryAlgorithm.snapshotTime = atof(lineVector[1].c_str());
            }
        },
        [] (const SimulConfig& sc) {
            return vector<string> { to_string(sc.chemParams.chemistryAlgorithm.snapshotTime) };
        }
    );
    chemParser.addStringArgsWithAliases(
        "MINIMIZATIONTIME", { "MINIMIZATIONTIME:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            if(lineVector.size() != 2) {
                cout <<
                     "There was an error parsing input file at Chemistry algorithm. Exiting."
                     << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                sc.chemParams.chemistryAlgorithm.minimizationTime = atof(lineVector[1].c_str());
            }
        },
        [] (const SimulConfig& sc) {
            return vector<string> { to_string(sc.chemParams.chemistryAlgorithm.minimizationTime) };
        }
    );
    chemParser.addStringArgsWithAliases(
        "NEIGHBORLISTTIME", { "NEIGHBORLISTTIME:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            if(lineVector.size() != 2) {
                cout <<
                     "There was an error parsing input file at Chemistry algorithm. Exiting."
                     << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                sc.chemParams.chemistryAlgorithm.neighborListTime = atof(lineVector[1].c_str());
            }
        },
        [] (const SimulConfig& sc) {
            return vector<string> { to_string(sc.chemParams.chemistryAlgorithm.neighborListTime) };
        }
    );
    chemParser.addStringArgsWithAliases(
        "DATADUMPTIME", { "DATADUMPTIME:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            if(lineVector.size() > 2) {
                cout <<
                     "There was an error parsing input file at Chemistry algorithm. Exiting."
                     << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                sc.chemParams.chemistryAlgorithm.datadumpTime = atof(lineVector[1].c_str());
            }
        },
        [] (const SimulConfig& sc) {
            return vector<string> { to_string(sc.chemParams.chemistryAlgorithm.datadumpTime) };
        }
    );
    chemParser.addEmptyLine();
    chemParser.addStringArgsWithAliases(
        "RUNSTEPS", { "RUNSTEPS:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            if(lineVector.size() > 2) {
                cout <<
                     "There was an error parsing input file at Chemistry algorithm. Exiting."
                     << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                sc.chemParams.chemistryAlgorithm.runSteps = atoi(lineVector[1].c_str());
            }
        },
        [] (const SimulConfig& sc) {
            return vector<string> { to_string(sc.chemParams.chemistryAlgorithm.runSteps) };
        }
    );
    chemParser.addStringArgsWithAliases(
        "SNAPSHOTSTEPS", { "SNAPSHOTSTEPS:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            if(lineVector.size() > 2) {
                cout <<
                     "There was an error parsing input file at Chemistry algorithm. Exiting."
                     << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                sc.chemParams.chemistryAlgorithm.snapshotSteps = atoi(lineVector[1].c_str());
            }
        },
        [] (const SimulConfig& sc) {
            return vector<string> { to_string(sc.chemParams.chemistryAlgorithm.snapshotSteps) };
        }
    );
    chemParser.addStringArgsWithAliases(
        "MINIMIZATIONSTEPS", { "MINIMIZATIONSTEPS:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            if(lineVector.size() != 2) {
                cout <<
                     "There was an error parsing input file at Chemistry algorithm. Exiting."
                     << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                sc.chemParams.chemistryAlgorithm.minimizationSteps = atoi(lineVector[1].c_str());
            }
        },
        [] (const SimulConfig& sc) {
            return vector<string> { to_string(sc.chemParams.chemistryAlgorithm.minimizationSteps) };
        }
    );
    chemParser.addStringArgsWithAliases(
        "NEIGHBORLISTSTEPS", { "NEIGHBORLISTSTEPS:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            if(lineVector.size() != 2) {
                cout <<
                     "There was an error parsing input file at Chemistry algorithm. Exiting."
                     << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                sc.chemParams.chemistryAlgorithm.neighborListSteps = atoi(lineVector[1].c_str());
            }
        },
        [] (const SimulConfig& sc) {
            return vector<string> { to_string(sc.chemParams.chemistryAlgorithm.neighborListSteps) };
        }
    );
    chemParser.addEmptyLine();

    chemParser.addComment("###### Chemistry setup ######");
    chemParser.addEmptyLine();
    chemParser.addStringArgsWithAliases(
        "CHEMISTRYFILE", { "CHEMISTRYFILE:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            if(lineVector.size() > 2) {
                cout << "Error reading chemistry input file. Exiting." << endl;
                exit(EXIT_FAILURE);
            }

            else if (lineVector.size() == 2)
                sc.chemParams.chemistrySetup.inputFile = lineVector[1];

            else if(lineVector.size() < 2) {
                cout << "Must specify a chemistry input file. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
        },
        [] (const SimulConfig& sc) {
            return vector<string> { sc.chemParams.chemistrySetup.inputFile.string() };
        }
    );
    chemParser.addEmptyLine();

    chemParser.addComment("###### Numbers and protocols ######");
    chemParser.addEmptyLine();

    chemParser.addComment("### Numbers");
    chemParser.addStringArgsWithAliases(
        "NUMFILAMENTTYPES", { "NUMFILAMENTTYPES:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            if(lineVector.size() > 2) {
                cout <<
                     "There was an error parsing input file at Chemistry parameters. Exiting."
                     << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                sc.chemParams.numFilaments = atoi(lineVector[1].c_str());
            }
        },
        [] (const SimulConfig& sc) {
            return vector<string> { to_string(sc.chemParams.numFilaments) };
        }
    );
    chemParser.addStringArgsWithAliases(
        "NUMBINDINGSITES", { "NUMBINDINGSITES:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            if (lineVector.size() >= 2) {
                for(int i = 1; i < lineVector.size(); i++)
                    sc.chemParams.numBindingSites.push_back(atoi(lineVector[i].c_str()));
            }
        },
        [] (const SimulConfig& sc) {
            vector<string> res;
            for(auto x : sc.chemParams.numBindingSites) {
                res.push_back(to_string(x));
            }
            return res;
        }
    );
    chemParser.addStringArgsWithAliases(
        "NUMMOTORHEADSMIN", { "NUMMOTORHEADSMIN:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            if (lineVector.size() >= 2) {
                for(int i = 1; i < lineVector.size(); i++)
                    sc.chemParams.motorNumHeadsMin.push_back(atoi(lineVector[i].c_str()));
            }
        },
        [] (const SimulConfig& sc) {
            vector<string> res;
            for(auto x : sc.chemParams.motorNumHeadsMin) {
                res.push_back(to_string(x));
            }
            return res;
        }
    );
    chemParser.addStringArgsWithAliases(
        "NUMMOTORHEADSMAX", { "NUMMOTORHEADSMAX:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            if (lineVector.size() >= 2) {
                for(int i = 1; i < lineVector.size(); i++)
                    sc.chemParams.motorNumHeadsMax.push_back(atoi(lineVector[i].c_str()));
            }
        },
        [] (const SimulConfig& sc) {
            vector<string> res;
            for(auto x : sc.chemParams.motorNumHeadsMax) {
                res.push_back(to_string(x));
            }
            return res;
        }
    );
    chemParser.addStringArgsWithAliases(
        "MOTORSTEPSIZE", { "MOTORSTEPSIZE:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            if (lineVector.size() >= 2) {
                for(int i = 1; i < lineVector.size(); i++)
                    sc.chemParams.motorStepSize.push_back(atof(lineVector[i].c_str()));
            }
        },
        [] (const SimulConfig& sc) {
            vector<string> res;
            for(auto x : sc.chemParams.motorStepSize) {
                res.push_back(to_string(x));
            }
            return res;
        }
    );
    chemParser.addStringArgsWithAliases(
        "LINKERBINDINGSKIP", { "LINKERBINDINGSKIP:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            if(lineVector.size() != 2) {
                cout <<
                "There was an error parsing input file at Chemistry algorithm. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                sc.chemParams.linkerbindingskip = atoi(lineVector[1].c_str());
            }
        },
        [] (const SimulConfig& sc) {
            return vector<string> { to_string(sc.chemParams.linkerbindingskip) };
        }
    );
    chemParser.addEmptyLine();

    chemParser.addComment("### Special protocols");
    chemParser.addStringArgsWithAliases(
        "SPECIALPROTOCOL", { "SPECIALPROTOCOL:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            //the vector size can be 5 for PINLOWERBOUNDARYFILAMENTS
            if(lineVector.size() > 7) {
                cout <<
                     "There was an error parsing input file at Chemistry parameters. Exiting."
                     << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 3) {
                if(lineVector[1] == "MAKELINKERSSTATIC") {
                    sc.chemParams.makeLinkersStatic = true;
                    sc.chemParams.makeLinkersStaticTime = atof(lineVector[2].c_str());
                }
                if(lineVector[1] == "MAKEFILAMENTSSTATIC") {
                    sc.chemParams.makeFilamentsStatic = true;
                    sc.chemParams.makeFilamentsStaticTime = atof(lineVector[2].c_str());
                }
                
            }
            else if (lineVector.size() == 4) {
                if(lineVector[1]  == "RATEDEPEND") {
                    sc.chemParams.makeRateDepend = true;
                    sc.chemParams.makeRateDependTime = atof(lineVector[2].c_str());
                    sc.chemParams.makeRateDependForce = atof(lineVector[3].c_str());
                }
            }
            else if (lineVector.size() == 7) {
                if(lineVector[1]  == "AFM") {
                    sc.chemParams.makeAFM = true;
                    //displacement of each pull
                    sc.chemParams.AFMStep1 = atof(lineVector[2].c_str());
                    sc.chemParams.AFMStep2 = atof(lineVector[3].c_str());
                    //change dispalcement from 1 to 2
                    sc.chemParams.IterChange = atof(lineVector[4].c_str());
                    //total step of each AFM pull
                    sc.chemParams.StepTotal = atof(lineVector[5].c_str());
                    //time between each pull
                    sc.chemParams.StepTime = atof(lineVector[6].c_str());
                }
            }
        },
        [] (const SimulConfig& sc) {
            vector<vector<string>> res;
            if(sc.chemParams.makeLinkersStatic) {
                res.push_back({ "MAKELINKERSSTATIC", to_string(sc.chemParams.makeLinkersStaticTime) });
            }
            if(sc.chemParams.makeFilamentsStatic) {
                res.push_back({ "MAKEFILAMENTSSTATIC", to_string(sc.chemParams.makeFilamentsStaticTime) });
            }
            if(sc.chemParams.makeRateDepend) {
                res.push_back({
                    "RATEDEPEND",
                    to_string(sc.chemParams.makeRateDependTime),
                    to_string(sc.chemParams.makeRateDependForce)
                });
            }
            if(sc.chemParams.makeAFM) {
                res.push_back({
                    "AFM",
                    to_string(sc.chemParams.AFMStep1),
                    to_string(sc.chemParams.AFMStep2),
                    to_string(sc.chemParams.IterChange),
                    to_string(sc.chemParams.StepTotal),
                    to_string(sc.chemParams.StepTime)
                });
            }
            return res;
        }
    );
    chemParser.addEmptyLine();

    chemParser.addComment("### Dissipation tracking");
    chemParser.addStringArgsWithAliases(
        "DISSIPATIONTRACKING", { "DISSIPATIONTRACKING:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            if(lineVector.size() != 2) {
                cout <<
                "There was an error parsing input file at Chemistry algorithm. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                sc.chemParams.dissTracking = (lineVector[1] == "ON");
            }
        },
        [] (const SimulConfig& sc) {
            return vector<string> { sc.chemParams.dissTracking ? "ON" : "OFF" };
        }
    );
    chemParser.addStringArgsWithAliases(
        "EVENTTRACKING", { "EVENTTRACKING:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            if(lineVector.size() != 2) {
                cout <<
                "There was an error parsing input file at Chemistry algorithm. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                sc.chemParams.eventTracking = (lineVector[1] == "ON");
            }
        },
        [] (const SimulConfig& sc) {
            return vector<string> { sc.chemParams.eventTracking ? "ON" : "OFF" };
        }
    );
    chemParser.addEmptyLine();

}

void SystemParser::chemPostProcessing(SimulConfig& sc) const {
    sc.chemParams.bindingSites.clear();

    //Figure out the binding sites
    for(int i = 0; i < sc.chemParams.numBindingSites.size(); i++) {

        sc.chemParams.maxbindingsitespercylinder = max(sc.chemParams.maxbindingsitespercylinder,
                                                 sc.chemParams.numBindingSites[i]);

        vector<short> tempBindingSites;

        int deltaBinding = sc.geoParams.cylinderNumMon[i] /
                           sc.chemParams.numBindingSites[i];

        int firstBindingSite = deltaBinding / 2 + 1;
        int bindingCount = firstBindingSite;

        //add all other binding sites
        while(bindingCount < sc.geoParams.cylinderNumMon[i]) {
            tempBindingSites.push_back(bindingCount);
            bindingCount += deltaBinding;
        }


        //push to sc.chemParams
        sc.chemParams.bindingSites.push_back(tempBindingSites);
    }

    //Find the maximum allowed Cindex and shift operator
    auto np2 = mathfunc::nextPowerOf2(uint32_t(sc.chemParams
            .maxbindingsitespercylinder));

    if(np2 == sc.chemParams.maxbindingsitespercylinder)
        np2 *= 2;

	sc.chemParams.shiftbybits = log2(np2);
    sc.chemParams.maxStableIndex = numeric_limits<uint32_t>::max()/sc.chemParams.shiftbybits -1;

}

void SystemParser::initMechParser() {
    using namespace std;

    mechParser.addComment("##################################################");
    mechParser.addComment("### Mechanical parameters");
    mechParser.addComment("##################################################");
    mechParser.addEmptyLine();

    mechParser.addComment("###### Mechanical algorithm ######");
    mechParser.addEmptyLine();

    mechParser.addStringArgsWithAliases(
        "CONJUGATEGRADIENT", { "CONJUGATEGRADIENT:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            if(lineVector.size() != 2) {
                cout <<
                     "A conjugate gradient method must be specified. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                sc.mechParams.mechanicsAlgorithm.ConjugateGradient = lineVector[1];
            }
        },
        [] (const SimulConfig& sc) {
            vector<vector<string>> res;
            if(const auto& s = sc.mechParams.mechanicsAlgorithm.ConjugateGradient; !s.empty()) {
                res.push_back({ s });
            }
            return res;
        }
    );
    mechParser.addStringArgsWithAliases(
        "MD", {},
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            if(lineVector.size() != 2) {
                cout << "A Mechanics algorithm must be specified. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                sc.mechParams.mechanicsAlgorithm.MD = lineVector[1];
            }
        },
        [] (const SimulConfig& sc) {
            vector<vector<string>> res;
            if(const auto& s = sc.mechParams.mechanicsAlgorithm.MD; !s.empty()) {
                res.push_back({ s });
            }
            return res;
        }
    );
    mechParser.addStringArgsWithAliases(
        "GRADIENTTOLERANCE", { "GRADIENTTOLERANCE:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            if (lineVector.size() == 2) {
                sc.mechParams.mechanicsAlgorithm.gradientTolerance = atof(lineVector[1].c_str());
            }
        },
        [] (const SimulConfig& sc) {
            vector<string> res;
            if(const auto& x = sc.mechParams.mechanicsAlgorithm.gradientTolerance; x) {
                res.push_back(to_string(x));
            }
            return res;
        }
    );
    mechParser.addStringArgsWithAliases(
        "MAXDISTANCE", { "MAXDISTANCE:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            if (lineVector.size() == 2) {
                sc.mechParams.mechanicsAlgorithm.maxDistance = atof(lineVector[1].c_str());
            }
        },
        [] (const SimulConfig& sc) {
            vector<string> res;
            if(const auto& x = sc.mechParams.mechanicsAlgorithm.maxDistance; x) {
                res.push_back(to_string(x));
            }
            return res;
        }
    );
    mechParser.addStringArgsWithAliases(
        "LAMBDAMAX", { "LAMBDAMAX:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            if (lineVector.size() == 2) {
                sc.mechParams.mechanicsAlgorithm.lambdaMax = atof(lineVector[1].c_str());
            }
        },
        [] (const SimulConfig& sc) {
            vector<string> res;
            if(const auto& x = sc.mechParams.mechanicsAlgorithm.lambdaMax; x) {
                res.push_back(to_string(x));
            }
            return res;
        }
    );
    mechParser.addStringArgsWithAliases(
        "LAMBDARUNNINGAVERAGEPROBABILITY", { "LAMBDARUNNINGAVERAGEPROBABILITY:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            if (lineVector.size() == 2) {
                sc.mechParams.mechanicsAlgorithm.lambdarunningaverageprobability = atof(lineVector[1].c_str());
            }
        },
        [] (const SimulConfig& sc) {
            vector<string> res;
            if(const auto& x = sc.mechParams.mechanicsAlgorithm.lambdarunningaverageprobability; x) {
                res.push_back(to_string(x));
            }
            return res;
        }
    );
    mechParser.addStringArgsWithAliases(
        "LINESEARCHALGORITHM", { "LINESEARCHALGORITHM:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            if (lineVector.size() == 2) {
                sc.mechParams.mechanicsAlgorithm.linesearchalgorithm = lineVector[1];
            }
        },
        [] (const SimulConfig& sc) {
            vector<vector<string>> res;
            if(const auto& s = sc.mechParams.mechanicsAlgorithm.linesearchalgorithm; !s.empty()) {
                res.push_back({ s });
            }
            return res;
        }
    );
    mechParser.addEmptyLine();

    mechParser.addComment("###### Force fields ######");
    mechParser.addEmptyLine();

    mechParser.addComment("### Actin filaments ");
    mechParser.addEmptyLine();

    mechParser.addComment("# Stretching: Popov et al, 2016, PLoS Comp Biol");
    mechParser.addStringArgsWithAliases(
        "FSTRETCHINGFFTYPE", { "FSTRETCHINGFFTYPE:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            if(lineVector.size() > 2) {
                cout <<
                     "There was an error parsing input file at Filament stretching FF type. Exiting."
                     << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                sc.mechParams.mechanicsFFType.FStretchingType = lineVector[1];
            }
        },
        [] (const SimulConfig& sc) {
            vector<string> res;
            if(const auto& s = sc.mechParams.mechanicsFFType.FStretchingType; !s.empty()) {
                res.push_back(s);
            }
            return res;
        }
    );
    mechParser.addStringArgsWithAliases(
        "FSTRETCHINGK", { "FSTRETCHINGK:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            sc.mechParams.FStretchingK.clear();
            if (lineVector.size() >= 2) {
                for(int i = 1; i < lineVector.size(); i++)
                    sc.mechParams.FStretchingK.push_back(atof((lineVector[i].c_str())));
            }
        },
        [] (const SimulConfig& sc) {
            vector<string> res;
            for(const auto& k : sc.mechParams.FStretchingK) {
                res.push_back(to_string(k));
            }
            return res;
        }
    );
    mechParser.addEmptyLine();

    mechParser.addComment("# Bending: Ott et al, 1993, Phys Rev E");
    mechParser.addStringArgsWithAliases(
        "FBENDINGFFTYPE", { "FBENDINGFFTYPE:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            if(lineVector.size() > 2) {
                cout <<
                     "There was an error parsing input file at Filament bending FF type. Exiting."
                     << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                sc.mechParams.mechanicsFFType.FBendingType = lineVector[1];
            }
        },
        [] (const SimulConfig& sc) {
            vector<string> res;
            if(const auto& s = sc.mechParams.mechanicsFFType.FBendingType; !s.empty()) {
                res.push_back(s);
            }
            return res;
        }
    );
    mechParser.addStringArgsWithAliases(
        "FBENDINGK", { "FBENDINGK:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            sc.mechParams.FBendingK.clear();
            if (lineVector.size() >= 2) {
                for(int i = 1; i < lineVector.size(); i++)
                    sc.mechParams.FBendingK.push_back(atof((lineVector[i].c_str())));
            }
        },
        [] (const SimulConfig& sc) {
            vector<string> res;
            for(const auto& k : sc.mechParams.FBendingK) {
                res.push_back(to_string(k));
            }
            return res;
        }
    );
    mechParser.addStringArgsWithAliases(
        "FBENDINGTHETA", { "FBENDINGTHETA:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            sc.mechParams.FBendingTheta.clear();
            if (lineVector.size() >= 2) {
                for(int i = 1; i < lineVector.size(); i++)
                    sc.mechParams.FBendingTheta.push_back(atof((lineVector[i].c_str())));
            }
        },
        [] (const SimulConfig& sc) {
            vector<string> res;
            for(const auto& k : sc.mechParams.FBendingTheta) {
                res.push_back(to_string(k));
            }
            return res;
        }
    );
    mechParser.addEmptyLine();

    mechParser.addComment("# Twisting: Currently not implemented");
    mechParser.addStringArgsWithAliases(
        "FTWISTINGFFTYPE", { "FTWISTINGFFTYPE:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            if(lineVector.size() > 2) {
                cout <<
                     "There was an error parsing input file at Filament twisting FF type. Exiting."
                     << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                sc.mechParams.mechanicsFFType.FTwistingType = lineVector[1];
            }
        },
        [] (const SimulConfig& sc) {
            vector<string> res;
            if(const auto& s = sc.mechParams.mechanicsFFType.FTwistingType; !s.empty()) {
                res.push_back(s);
            }
            return res;
        }
    );
    mechParser.addStringArgsWithAliases(
        "FTWISTINGK", { "FTWISTINGK:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            sc.mechParams.FTwistingK.clear();
            if (lineVector.size() >= 2) {
                for(int i = 1; i < lineVector.size(); i++)
                    sc.mechParams.FTwistingK.push_back(atof((lineVector[i].c_str())));
            }
        },
        [] (const SimulConfig& sc) {
            vector<string> res;
            for(const auto& k : sc.mechParams.FTwistingK) {
                res.push_back(to_string(k));
            }
            return res;
        }
    );
    mechParser.addStringArgsWithAliases(
        "FTWISTINGPHI", { "FTWISTINGPHI:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            sc.mechParams.FTwistingPhi.clear();
            if (lineVector.size() >= 2) {
                for(int i = 1; i < lineVector.size(); i++)
                    sc.mechParams.FTwistingPhi.push_back(atof((lineVector[i].c_str())));
            }
        },
        [] (const SimulConfig& sc) {
            vector<string> res;
            for(const auto& k : sc.mechParams.FTwistingPhi) {
                res.push_back(to_string(k));
            }
            return res;
        }
    );
    mechParser.addEmptyLine();

    mechParser.addComment("### Linkers");
    mechParser.addEmptyLine();
    mechParser.addComment("# Stretching: Didonna et al, Phys Rev E, 2007");
    mechParser.addStringArgsWithAliases(
        "LSTRETCHINGFFTYPE", { "LSTRETCHINGFFTYPE:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            if(lineVector.size() > 2) {
                cout <<
                     "There was an error parsing input file at Linker stretching FF type. Exiting."
                     << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                sc.mechParams.mechanicsFFType.LStretchingType = lineVector[1];
            }
        },
        [] (const SimulConfig& sc) {
            vector<string> res;
            if(const auto& s = sc.mechParams.mechanicsFFType.LStretchingType; !s.empty()) {
                res.push_back(s);
            }
            return res;
        }
    );
    mechParser.addStringArgsWithAliases(
        "LSTRETCHINGK", { "LSTRETCHINGK:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            sc.mechParams.LStretchingK.clear();
            if (lineVector.size() >= 2) {
                for(int i = 1; i < lineVector.size(); i++)
                    sc.mechParams.LStretchingK.push_back(atof((lineVector[i].c_str())));
            }
        },
        [] (const SimulConfig& sc) {
            vector<string> res;
            for(const auto& k : sc.mechParams.LStretchingK) {
                res.push_back(to_string(k));
            }
            return res;
        }
    );
    mechParser.addEmptyLine();
    mechParser.addComment("# Bending/Twisting: not implemented");
    mechParser.addStringArgsWithAliases(
        "LBENDINGFFTYPE", { "LBENDINGFFTYPE:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            if(lineVector.size() > 2) {
                cout <<
                     "There was an error parsing input file at Linker bending FF type. Exiting."
                     << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                sc.mechParams.mechanicsFFType.LBendingType = lineVector[1];
            }
        },
        [] (const SimulConfig& sc) {
            vector<string> res;
            if(const auto& s = sc.mechParams.mechanicsFFType.LBendingType; !s.empty()) {
                res.push_back(s);
            }
            return res;
        }
    );
    mechParser.addStringArgsWithAliases(
        "LBENDINGK", { "LBENDINGK:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            sc.mechParams.LBendingK.clear();
            if (lineVector.size() >= 2) {
                for(int i = 1; i < lineVector.size(); i++)
                    sc.mechParams.LBendingK.push_back(atof((lineVector[i].c_str())));
            }
        },
        [] (const SimulConfig& sc) {
            vector<string> res;
            for(const auto& k : sc.mechParams.LBendingK) {
                res.push_back(to_string(k));
            }
            return res;
        }
    );
    mechParser.addStringArgsWithAliases(
        "LBENDINGTHETA", { "LBENDINGTHETA:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            sc.mechParams.LBendingTheta.clear();
            if (lineVector.size() >= 2) {
                for(int i = 1; i < lineVector.size(); i++)
                    sc.mechParams.LBendingTheta.push_back(atof((lineVector[i].c_str())));
            }
        },
        [] (const SimulConfig& sc) {
            vector<string> res;
            for(const auto& k : sc.mechParams.LBendingTheta) {
                res.push_back(to_string(k));
            }
            return res;
        }
    );
    mechParser.addStringArgsWithAliases(
        "LTWISTINGFFTYPE", { "LTWISTINGFFTYPE:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            if(lineVector.size() > 2) {
                cout <<
                     "There was an error parsing input file at Linker twisting FF type. Exiting."
                     << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                sc.mechParams.mechanicsFFType.LTwistingType = lineVector[1];
            }
        },
        [] (const SimulConfig& sc) {
            vector<string> res;
            if(const auto& s = sc.mechParams.mechanicsFFType.LTwistingType; !s.empty()) {
                res.push_back(s);
            }
            return res;
        }
    );
    mechParser.addStringArgsWithAliases(
        "LTWISTINGK", { "LTWISTINGK:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            sc.mechParams.LTwistingK.clear();
            if (lineVector.size() >= 2) {
                for(int i = 1; i < lineVector.size(); i++)
                    sc.mechParams.LTwistingK.push_back(atof((lineVector[i].c_str())));
            }
        },
        [] (const SimulConfig& sc) {
            vector<string> res;
            for(const auto& k : sc.mechParams.LTwistingK) {
                res.push_back(to_string(k));
            }
            return res;
        }
    );
    mechParser.addStringArgsWithAliases(
        "LTWISTINGPHI", { "LTWISTINGPHI:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            sc.mechParams.LTwistingPhi.clear();
            if (lineVector.size() >= 2) {
                for(int i = 1; i < lineVector.size(); i++)
                    sc.mechParams.LTwistingPhi.push_back(atof((lineVector[i].c_str())));
            }
        },
        [] (const SimulConfig& sc) {
            vector<string> res;
            for(const auto& k : sc.mechParams.LTwistingPhi) {
                res.push_back(to_string(k));
            }
            return res;
        }
    );
    mechParser.addEmptyLine();

    mechParser.addComment("### Myosin motors");
    mechParser.addEmptyLine();
    mechParser.addComment("# Stretching: Vilfan, Biophys J, 2010");
    mechParser.addStringArgsWithAliases(
        "MSTRETCHINGFFTYPE", { "MSTRETCHINGFFTYPE:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            if(lineVector.size() > 2) {
                cout <<
                     "There was an error parsing input file at Motor stretching FF type. Exiting."
                     << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                sc.mechParams.mechanicsFFType.MStretchingType = lineVector[1];
            }
        },
        [] (const SimulConfig& sc) {
            vector<string> res;
            if(const auto& s = sc.mechParams.mechanicsFFType.MStretchingType; !s.empty()) {
                res.push_back(s);
            }
            return res;
        }
    );
    mechParser.addStringArgsWithAliases(
        "MSTRETCHINGK", { "MSTRETCHINGK:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            sc.mechParams.MStretchingK.clear();
            if (lineVector.size() >= 2) {
                for(int i = 1; i < lineVector.size(); i++)
                    sc.mechParams.MStretchingK.push_back(atof((lineVector[i].c_str())));
            }
        },
        [] (const SimulConfig& sc) {
            vector<string> res;
            for(const auto& k : sc.mechParams.MStretchingK) {
                res.push_back(to_string(k));
            }
            return res;
        }
    );
    mechParser.addEmptyLine();
    mechParser.addComment("# Bending/Twisting: not implemented");
    mechParser.addStringArgsWithAliases(
        "MBENDINGFFTYPE", { "MBENDINGFFTYPE:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            if(lineVector.size() > 2) {
                cout <<
                     "There was an error parsing input file at Motor bending FF type. Exiting."
                     << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                sc.mechParams.mechanicsFFType.MBendingType = lineVector[1];
            }
        },
        [] (const SimulConfig& sc) {
            vector<string> res;
            if(const auto& s = sc.mechParams.mechanicsFFType.MBendingType; !s.empty()) {
                res.push_back(s);
            }
            return res;
        }
    );
    mechParser.addStringArgsWithAliases(
        "MBENDINGK", { "MBENDINGK:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            sc.mechParams.MBendingK.clear();
            if (lineVector.size() >= 2) {
                for(int i = 1; i < lineVector.size(); i++)
                    sc.mechParams.MBendingK.push_back(atof((lineVector[i].c_str())));
            }
        },
        [] (const SimulConfig& sc) {
            vector<string> res;
            for(const auto& k : sc.mechParams.MBendingK) {
                res.push_back(to_string(k));
            }
            return res;
        }
    );
    mechParser.addStringArgsWithAliases(
        "MBENDINGTHETA", { "MBENDINGTHETA:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            sc.mechParams.MBendingTheta.clear();
            if (lineVector.size() >= 2) {
                for(int i = 1; i < lineVector.size(); i++)
                    sc.mechParams.MBendingTheta.push_back(atof((lineVector[i].c_str())));
            }
        },
        [] (const SimulConfig& sc) {
            vector<string> res;
            for(const auto& k : sc.mechParams.MBendingTheta) {
                res.push_back(to_string(k));
            }
            return res;
        }
    );
    mechParser.addStringArgsWithAliases(
        "MTWISTINGFFTYPE", { "MTWISTINGFFTYPE:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            if(lineVector.size() > 2) {
                cout <<
                     "There was an error parsing input file at Motor twisting FF type. Exiting."
                     << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                sc.mechParams.mechanicsFFType.MTwistingType = lineVector[1];
            }
        },
        [] (const SimulConfig& sc) {
            vector<string> res;
            if(const auto& s = sc.mechParams.mechanicsFFType.MTwistingType; !s.empty()) {
                res.push_back(s);
            }
            return res;
        }
    );
    mechParser.addStringArgsWithAliases(
        "MTWISTINGK", { "MTWISTINGK:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            sc.mechParams.MTwistingK.clear();
            if (lineVector.size() >= 2) {
                for(int i = 1; i < lineVector.size(); i++)
                    sc.mechParams.MTwistingK.push_back(atof((lineVector[i].c_str())));
            }
        },
        [] (const SimulConfig& sc) {
            vector<string> res;
            for(const auto& k : sc.mechParams.MTwistingK) {
                res.push_back(to_string(k));
            }
            return res;
        }
    );
    mechParser.addStringArgsWithAliases(
        "MTWISTINGPHI", { "MTWISTINGPHI:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            sc.mechParams.MTwistingPhi.clear();
            if (lineVector.size() >= 2) {
                for(int i = 1; i < lineVector.size(); i++)
                    sc.mechParams.MTwistingPhi.push_back(atof((lineVector[i].c_str())));
            }
        },
        [] (const SimulConfig& sc) {
            vector<string> res;
            for(const auto& k : sc.mechParams.MTwistingPhi) {
                res.push_back(to_string(k));
            }
            return res;
        }
    );
    mechParser.addEmptyLine();

    mechParser.addComment("### Arp2/3 brancher");
    mechParser.addComment("# 4 force fields: Popov et al, Plos Comp Biol, 2016");
    mechParser.addComment("# No reliable literature values for this FF");
    mechParser.addEmptyLine();
    mechParser.addStringArgsWithAliases(
        "BRSTRETCHINGFFTYPE", { "BRSTRETCHINGFFTYPE:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            if(lineVector.size() > 2) {
                cout <<
                     "There was an error parsing input file at Brancher stretching FF type. Exiting."
                     << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                sc.mechParams.mechanicsFFType.BrStretchingType = lineVector[1];
            }
        },
        [] (const SimulConfig& sc) {
            vector<string> res;
            if(const auto& s = sc.mechParams.mechanicsFFType.BrStretchingType; !s.empty()) {
                res.push_back(s);
            }
            return res;
        }
    );
    mechParser.addStringArgsWithAliases(
        "BRSTRETCHINGK", { "BRSTRETCHINGK:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            sc.mechParams.BrStretchingK.clear();
            if (lineVector.size() >= 2) {
                for(int i = 1; i < lineVector.size(); i++)
                    sc.mechParams.BrStretchingK.push_back(atof((lineVector[i].c_str())));
            }
        },
        [] (const SimulConfig& sc) {
            vector<string> res;
            for(const auto& k : sc.mechParams.BrStretchingK) {
                res.push_back(to_string(k));
            }
            return res;
        }
    );
    mechParser.addStringArgsWithAliases(
        "BRSTRETCHINGL", { "BRSTRETCHINGL:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            sc.mechParams.BrStretchingL.clear();
            if (lineVector.size() >= 2) {
                for(int i = 1; i < lineVector.size(); i++)
                    sc.mechParams.BrStretchingL.push_back(atof((lineVector[i].c_str())));
            }
        },
        [] (const SimulConfig& sc) {
            vector<string> res;
            for(const auto& k : sc.mechParams.BrStretchingL) {
                res.push_back(to_string(k));
            }
            return res;
        }
    );
    mechParser.addEmptyLine();
    mechParser.addStringArgsWithAliases(
        "BRBENDINGFFTYPE", { "BRBENDINGFFTYPE:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            if(lineVector.size() > 2) {
                cout <<
                     "There was an error parsing input file at Brancher bending FF type. Exiting."
                     << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                sc.mechParams.mechanicsFFType.BrBendingType = lineVector[1];
            }
        },
        [] (const SimulConfig& sc) {
            vector<string> res;
            if(const auto& s = sc.mechParams.mechanicsFFType.BrBendingType; !s.empty()) {
                res.push_back(s);
            }
            return res;
        }
    );
    mechParser.addStringArgsWithAliases(
        "BRBENDINGK", { "BRBENDINGK:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            sc.mechParams.BrBendingK.clear();
            if (lineVector.size() >= 2) {
                for(int i = 1; i < lineVector.size(); i++)
                    sc.mechParams.BrBendingK.push_back(atof((lineVector[i].c_str())));
            }
        },
        [] (const SimulConfig& sc) {
            vector<string> res;
            for(const auto& k : sc.mechParams.BrBendingK) {
                res.push_back(to_string(k));
            }
            return res;
        }
    );
    mechParser.addStringArgsWithAliases(
        "BRBENDINGTHETA", { "BRBENDINGTHETA:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            sc.mechParams.BrBendingTheta.clear();
            if (lineVector.size() >= 2) {
                for(int i = 1; i < lineVector.size(); i++)
                    sc.mechParams.BrBendingTheta.push_back(atof((lineVector[i].c_str())));
            }
        },
        [] (const SimulConfig& sc) {
            vector<string> res;
            for(const auto& k : sc.mechParams.BrBendingTheta) {
                res.push_back(to_string(k));
            }
            return res;
        }
    );
    mechParser.addEmptyLine();
    mechParser.addStringArgsWithAliases(
        "BRDIHEDRALFFTYPE", { "BRDIHEDRALFFTYPE:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            if(lineVector.size() > 2) {
                cout <<
                     "There was an error parsing input file at Brancher dihedral FF type. Exiting."
                     << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                sc.mechParams.mechanicsFFType.BrDihedralType = lineVector[1];
            }
        },
        [] (const SimulConfig& sc) {
            vector<string> res;
            if(const auto& s = sc.mechParams.mechanicsFFType.BrDihedralType; !s.empty()) {
                res.push_back(s);
            }
            return res;
        }
    );
    mechParser.addStringArgsWithAliases(
        "BRDIHEDRALK", { "BRDIHEDRALK:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            sc.mechParams.BrDihedralK.clear();
            if (lineVector.size() >= 2) {
                for(int i = 1; i < lineVector.size(); i++)
                    sc.mechParams.BrDihedralK.push_back(atof((lineVector[i].c_str())));
            }
        },
        [] (const SimulConfig& sc) {
            vector<string> res;
            for(const auto& k : sc.mechParams.BrDihedralK) {
                res.push_back(to_string(k));
            }
            return res;
        }
    );
    mechParser.addEmptyLine();
    mechParser.addStringArgsWithAliases(
        "BRPOSITIONFFTYPE", { "BRPOSITIONFFTYPE:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            if(lineVector.size() > 2) {
                cout <<
                     "There was an error parsing input file at Brancher position FF type. Exiting."
                     << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                sc.mechParams.mechanicsFFType.BrPositionType = lineVector[1];
            }
        },
        [] (const SimulConfig& sc) {
            vector<string> res;
            if(const auto& s = sc.mechParams.mechanicsFFType.BrPositionType; !s.empty()) {
                res.push_back(s);
            }
            return res;
        }
    );
    mechParser.addStringArgsWithAliases(
        "BRPOSITIONK", { "BRPOSITIONK:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            sc.mechParams.BrPositionK.clear();
            if (lineVector.size() >= 2) {
                for(int i = 1; i < lineVector.size(); i++)
                    sc.mechParams.BrPositionK.push_back(atof((lineVector[i].c_str())));
            }
        },
        [] (const SimulConfig& sc) {
            vector<string> res;
            for(const auto& k : sc.mechParams.BrPositionK) {
                res.push_back(to_string(k));
            }
            return res;
        }
    );
    mechParser.addEmptyLine();

    mechParser.addComment("### Boundary");
    mechParser.addComment("# Boundary repulsion force field. Parameters are set in the boundary section.");
    mechParser.addStringArgsWithAliases(
        "BOUNDARYFFTYPE", { "BOUNDARYFFTYPE:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            if(lineVector.size() > 2) {
                cout <<
                     "There was an error parsing input file at Boundary FF type. Exiting."
                     << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                sc.mechParams.mechanicsFFType.BoundaryFFType = lineVector[1];
            }
        },
        [] (const SimulConfig& sc) {
            vector<string> res;
            if(const auto& s = sc.mechParams.mechanicsFFType.BoundaryFFType; !s.empty()) {
                res.push_back(s);
            }
            return res;
        }
    );
    mechParser.addEmptyLine();

    mechParser.addComment("### Volume exclusion: Popov et al, 2016, PLoS Comp Biol");
    mechParser.addEmptyLine();
    mechParser.addStringArgsWithAliases(
        "VOLUMEFFTYPE", { "VOLUMEFFTYPE:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            if(lineVector.size() > 2) {
                cout <<
                     "There was an error parsing input file at Volume FF type. Exiting."
                     << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                sc.mechParams.mechanicsFFType.VolumeFFType = lineVector[1];
            }
        },
        [] (const SimulConfig& sc) {
            vector<string> res;
            if(const auto& s = sc.mechParams.mechanicsFFType.VolumeFFType; !s.empty()) {
                res.push_back(s);
            }
            return res;
        }
    );
    mechParser.addStringArgsWithAliases(
        "VOLUMEK", { "VOLUMEK:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            sc.mechParams.VolumeK.clear();
            if (lineVector.size() >= 2) {
                for(int i = 1; i < lineVector.size(); i++)
                    sc.mechParams.VolumeK.push_back(atof((lineVector[i].c_str())));
            }
        },
        [] (const SimulConfig& sc) {
            vector<string> res;
            for(const auto& k : sc.mechParams.VolumeK) {
                res.push_back(to_string(k));
            }
            return res;
        }
    );
    mechParser.addStringArgsWithAliases(
        "VOLUMECUTOFF", { "VOLUMECUTOFF:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            if(lineVector.size() != 2) {
                cout <<
                     "Error reading Volume cutoff. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
            else {
                sc.mechParams.VolumeCutoff = stod(lineVector[1]);
            }
        },
        [] (const SimulConfig& sc) {
            vector<vector<string>> res;
            if(sc.mechParams.VolumeCutoff) {
                res.push_back({ to_string(sc.mechParams.VolumeCutoff) });
            }
            return res;
        }
    );
    mechParser.addEmptyLine();

    mechParser.addComment("### Bubble interactions");
    mechParser.addStringArgsWithAliases(
        "BUBBLEFFTYPE", { "BUBBLEFFTYPE:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            if(lineVector.size() > 2) {
                cout <<
                     "There was an error parsing input file at Bubble FF type. Exiting."
                     << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                sc.mechParams.mechanicsFFType.BubbleFFType = lineVector[1];
            }
        },
        [] (const SimulConfig& sc) {
            vector<string> res;
            if(const auto& s = sc.mechParams.mechanicsFFType.BubbleFFType; !s.empty()) {
                res.push_back(s);
            }
            return res;
        }
    );
    mechParser.addStringArgsWithAliases(
        "BUBBLEINTERACTIONK", { "BUBBLEINTERACTIONK:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            sc.mechParams.BubbleK.clear();
            if (lineVector.size() >= 2) {
                for(int i = 1; i < lineVector.size(); i++)
                    sc.mechParams.BubbleK.push_back(atof((lineVector[i].c_str())));
            }
        },
        [] (const SimulConfig& sc) {
            vector<string> res;
            for(const auto& k : sc.mechParams.BubbleK) {
                res.push_back(to_string(k));
            }
            return res;
        }
    );
    mechParser.addStringArgsWithAliases(
        "BUBBLESCREENLENGTH", { "BUBBLESCREENLENGTH:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            sc.mechParams.BubbleScreenLength.clear();
            if (lineVector.size() >= 2) {
                for(int i = 1; i < lineVector.size(); i++)
                    sc.mechParams.BubbleScreenLength.push_back(atof((lineVector[i].c_str())));
            }
        },
        [] (const SimulConfig& sc) {
            vector<string> res;
            for(const auto& k : sc.mechParams.BubbleScreenLength) {
                res.push_back(to_string(k));
            }
            return res;
        }
    );
    mechParser.addStringArgsWithAliases(
        "BUBBLECUTOFF", { "BUBBLECUTOFF:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            if(lineVector.size() != 2) {
                cout <<
                     "Error reading Bubble cutoff. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
            else {
                sc.mechParams.BubbleCutoff = stod(lineVector[1]);
            }
        },
        [] (const SimulConfig& sc) {
            vector<vector<string>> res;
            if(sc.mechParams.BubbleCutoff) {
                res.push_back({ to_string(sc.mechParams.BubbleCutoff) });
            }
            return res;
        }
    );
    mechParser.addStringArgsWithAliases(
        "BUBBLERADIUS", { "BUBBLERADIUS:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            sc.mechParams.BubbleRadius.clear();
            if (lineVector.size() >= 2) {
                for(int i = 1; i < lineVector.size(); i++)
                    sc.mechParams.BubbleRadius.push_back(atof((lineVector[i].c_str())));
            }
        },
        [] (const SimulConfig& sc) {
            vector<string> res;
            for(const auto& k : sc.mechParams.BubbleRadius) {
                res.push_back(to_string(k));
            }
            return res;
        }
    );
    mechParser.addStringArgsWithAliases(
        "NUMBUBBLETYPES", { "NUMBUBBLETYPES:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            if(lineVector.size() != 2) {
                cout <<
                     "Error reading number of Bubble types. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
            else {
                sc.mechParams.numBubbleTypes = atoi(lineVector[1].c_str());
            }
        },
        [] (const SimulConfig& sc) {
            vector<vector<string>> res;
            if(sc.mechParams.numBubbleTypes) {
                res.push_back({ to_string(sc.mechParams.numBubbleTypes) });
            }
            return res;
        }
    );
    mechParser.addEmptyLine();

    mechParser.addComment("### MTOC");
    mechParser.addStringArgsWithAliases(
        "MTOCFFTYPE", { "MTOCFFTYPE:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            if(lineVector.size() > 2) {
                cout <<
                     "There was an error parsing input file at MTOC FF type. Exiting."
                     << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                sc.mechParams.mechanicsFFType.MTOCFFType = lineVector[1];
            }
        },
        [] (const SimulConfig& sc) {
            vector<string> res;
            if(const auto& s = sc.mechParams.mechanicsFFType.MTOCFFType; !s.empty()) {
                res.push_back(s);
            }
            return res;
        }
    );
    mechParser.addStringArgsWithAliases(
        "MTOCBENDINGK", { "MTOCBENDINGK:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            sc.mechParams.MTOCBendingK.clear();
            if (lineVector.size() >= 2) {
                for(int i = 1; i < lineVector.size(); i++)
                    sc.mechParams.MTOCBendingK.push_back(atof((lineVector[i].c_str())));
            }
        },
        [] (const SimulConfig& sc) {
            vector<string> res;
            for(const auto& k : sc.mechParams.MTOCBendingK) {
                res.push_back(to_string(k));
            }
            return res;
        }
    );
    mechParser.addEmptyLine();

    mechParser.addComment("### AFM bead");
    mechParser.addStringArgsWithAliases(
        "AFMFFTYPE", { "AFMFFTYPE:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            if(lineVector.size() > 2) {
                cout <<
                     "There was an error parsing input file at AFM FF type. Exiting."
                     << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                sc.mechParams.mechanicsFFType.AFMFFType = lineVector[1];
            }
        },
        [] (const SimulConfig& sc) {
            vector<string> res;
            if(const auto& s = sc.mechParams.mechanicsFFType.AFMFFType; !s.empty()) {
                res.push_back(s);
            }
            return res;
        }
    );
    mechParser.addStringArgsWithAliases(
        "AFMBENDINGK", { "AFMBENDINGK:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            sc.mechParams.AFMBendingK.clear();
            if (lineVector.size() >= 2) {
                for(int i = 1; i < lineVector.size(); i++)
                    sc.mechParams.AFMBendingK.push_back(atof((lineVector[i].c_str())));
            }
        },
        [] (const SimulConfig& sc) {
            vector<string> res;
            for(const auto& k : sc.mechParams.AFMBendingK) {
                res.push_back(to_string(k));
            }
            return res;
        }
    );
    mechParser.addEmptyLine();

    mechParser.addComment("###### Protocols ######");
    mechParser.addEmptyLine();

    mechParser.addComment("# Hessian tracking");
    mechParser.addStringArgsWithAliases(
        "HESSIANTRACKING", { "HESSIANTRACKING:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            if(lineVector.size() != 3) {
                cout <<
                "There was an error parsing input file at Hessian tracking. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 3) {
                sc.mechParams.hessTracking = true;
                //sc.mechParams.hessDelta = atof(lineVector[1].c_str());
                sc.mechParams.hessSkip = atof(lineVector[1].c_str());
                int dense = atoi(lineVector[2].c_str());
                if(dense == 0){
                    sc.mechParams.denseEstimation = true;
                }else{
                    sc.mechParams.denseEstimation = false;
                }
                
            }
        },
        [] (const SimulConfig& sc) {
            vector<vector<string>> res;
            if(sc.mechParams.hessTracking) {
                res.push_back({
                    to_string(sc.mechParams.hessSkip),
                    to_string(sc.mechParams.denseEstimation ? 0 : 1)
                });
            }
            return res;
        }
    );
    mechParser.addEmptyLine();

    mechParser.addComment("# Same filament binding skip");
    mechParser.addStringArgsWithAliases(
        "SAMEFILBINDINGSKIP", { "SAMEFILBINDINGSKIP:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            if(lineVector.size() != 2) {
                cout <<
                "There was an error parsing input file at same filament binding skip. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                sc.mechParams.sameFilBindSkip = atoi(lineVector[1].c_str());
            }
        },
        [] (const SimulConfig& sc) {
            return vector<string> { to_string(sc.mechParams.sameFilBindSkip) };
        }
    );
    mechParser.addEmptyLine();

    mechParser.addComment("# Special protocols");
    mechParser.addStringArgsWithAliases(
        "SPECIALPROTOCOL", { "SPECIALPROTOCOL:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            if(lineVector[1] == "PINBOUNDARYFILAMENTS") {
                if(lineVector.size() > 5) {
                    cout <<
                         "There was an error parsing input file at pinning boundary filaments. Exiting."
                         << endl;
                    exit(EXIT_FAILURE);
                }

                else{
                    sc.mechParams.pinBoundaryFilaments = true;
                    sc.mechParams.pinK = atof(lineVector[2].c_str());
                    sc.mechParams.pinDistance = atof(lineVector[3].c_str());
                    sc.mechParams.pinTime = atof(lineVector[4].c_str());
                }

            }
            else if (lineVector.size() == 5) {

                //Qin
                if(lineVector[1] == "PINLOWERBOUNDARYFILAMENTS") {

                    sc.mechParams.pinLowerBoundaryFilaments = true;
                    sc.mechParams.pinK = atof(lineVector[2].c_str());
                    sc.mechParams.pinTime = atof(lineVector[3].c_str());
                    sc.mechParams.pinFraction = atof(lineVector[4].c_str());
                }
            }
        },
        [] (const SimulConfig& sc) {
            vector<vector<string>> res;
            if(sc.mechParams.pinBoundaryFilaments) {
                res.push_back({
                    "PINBOUNDARYFILAMENTS",
                    to_string(sc.mechParams.pinK),
                    to_string(sc.mechParams.pinDistance),
                    to_string(sc.mechParams.pinTime)
                });
            }
            if(sc.mechParams.pinLowerBoundaryFilaments) {
                res.push_back({
                    "PINLOWERBOUNDARYFILAMENTS",
                    to_string(sc.mechParams.pinK),
                    to_string(sc.mechParams.pinTime),
                    to_string(sc.mechParams.pinFraction)
                });
            }
            return res;
        }
    );
    mechParser.addEmptyLine();

}

void SystemParser::initBoundParser() {
    using namespace std;

    boundParser.addComment("##################################################");
    boundParser.addComment("### Boundary parameters");
    boundParser.addComment("##################################################");
    boundParser.addEmptyLine();

    boundParser.addComment("# Define network boundary geometry (CUBIC, SPHERICAL, CYLINDER)");
    boundParser.addEmptyLine();

    boundParser.addStringArgsWithAliases(
        "BOUNDARYSHAPE", { "BOUNDARYSHAPE:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            if(lineVector.size() != 2) {
                cout << "A boundary shape needs to be specified. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
            else {
                sc.boundParams.boundaryType.boundaryShape = lineVector[1];
            }
        },
        [] (const SimulConfig& sc) {
            return vector<string> { sc.boundParams.boundaryType.boundaryShape };
        }
    );
    boundParser.addEmptyLine();

    boundParser.addComment("# Define how network boundary moves");
    boundParser.addComment("# Usage: BOUNDARYMOVE LEFT/RIGHT/BOTTOM/TOP/FRONT/BACK");
    boundParser.addComment("# Changes are not recommended.");
    boundParser.addEmptyLine();

    boundParser.addStringArgsWithAliases(
        "BOUNDARYMOVE", { "BOUNDARYMOVE:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            if(lineVector.size() != 2) {
                cout << "A boundary move type needs to be specified. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
            else {
                sc.boundParams.boundaryType.boundaryMove.push_back(lineVector[1]);
            }
        },
        [] (const SimulConfig& sc) {
            vector<vector<string>> res;
            for(const auto& s : sc.boundParams.boundaryType.boundaryMove) {
                res.push_back({s});
            }
            return res;
        }
    );
    boundParser.addEmptyLine();

    boundParser.addComment("# Set diameter for SPHERICAL or CYLINDER type");
    boundParser.addComment("# CUBIC: No need to set, boundary is the same as network size");
    boundParser.addEmptyLine();

    boundParser.addStringArgsWithAliases(
        "BOUNDARYDIAMETER", { "BOUNDARYDIAMETER:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            if (lineVector.size() == 2) {
                sc.boundParams.diameter = atof((lineVector[1].c_str()));
            }
        },
        [] (const SimulConfig& sc) {
            return sc.boundParams.diameter ?
                vector<string> { to_string(sc.boundParams.diameter) } :
                vector<string> {};
        }
    );
    boundParser.addEmptyLine();

    boundParser.addComment("### Boundary interactions ");
    boundParser.addComment("# Repulsion: Popov et al, 2016, PLoS Comp Biol");
    boundParser.addEmptyLine();

    boundParser.addStringArgsWithAliases(
        "BOUNDARYCUTOFF", { "BOUNDARYCUTOFF:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            if(lineVector.size() != 2) {
                cout <<
                     "There was an error parsing input file at Boundary parameters. Exiting."
                     << endl;
                exit(EXIT_FAILURE);
            }
            else {
                sc.boundParams.BoundaryCutoff = atof((lineVector[1].c_str()));
            }
        },
        [] (const SimulConfig& sc) {
            return vector<string> { to_string(sc.boundParams.BoundaryCutoff) };
        }
    );
    boundParser.addStringArgsWithAliases(
        "BOUNDARYINTERACTIONK", { "BOUNDARYINTERACTIONK:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            if(lineVector.size() != 2) {
                cout <<
                     "There was an error parsing input file at Boundary parameters. Exiting."
                     << endl;
                exit(EXIT_FAILURE);
            }
            else {
                sc.boundParams.BoundaryK = atof((lineVector[1].c_str()));
            }
        },
        [] (const SimulConfig& sc) {
            return vector<string> { to_string(sc.boundParams.BoundaryK) };
        }
    );
    boundParser.addStringArgsWithAliases(
        "BOUNDARYSCREENLENGTH", { "BOUNDARYSCREENLENGTH:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            if(lineVector.size() != 2) {
                cout <<
                     "There was an error parsing input file at Boundary parameters. Exiting."
                     << endl;
                exit(EXIT_FAILURE);
            }
            else {
                sc.boundParams.BScreenLength = atof((lineVector[1].c_str()));
            }
        },
        [] (const SimulConfig& sc) {
            return vector<string> { to_string(sc.boundParams.BScreenLength) };
        }
    );
    boundParser.addEmptyLine();

    boundParser.addComment("# Set boundary move parameters");
    boundParser.addComment("# Changes are not recommended.");
    boundParser.addEmptyLine();

    boundParser.addStringArgsWithAliases(
        "BMOVESPEED", { "BMOVESPEED:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            if (lineVector.size() == 2) {
                sc.boundParams.moveSpeed = atof((lineVector[1].c_str()));
            }
        },
        [] (const SimulConfig& sc) {
            return sc.boundParams.moveSpeed ?
                vector<string> { to_string(sc.boundParams.moveSpeed) } :
                vector<string> {};
        }
    );
    boundParser.addStringArgsWithAliases(
        "BMOVESTARTTIME", { "BMOVESTARTTIME:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            if (lineVector.size() == 2) {
                sc.boundParams.moveStartTime = atof((lineVector[1].c_str()));
            }
        },
        [] (const SimulConfig& sc) {
            return sc.boundParams.moveStartTime ?
                vector<string> { to_string(sc.boundParams.moveStartTime) } :
                vector<string> {};
        }
    );
    boundParser.addStringArgsWithAliases(
        "BMOVEENDTIME", { "BMOVEENDTIME:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            if (lineVector.size() == 2) {
                sc.boundParams.moveEndTime = atof((lineVector[1].c_str()));
            }
        },
        [] (const SimulConfig& sc) {
            return sc.boundParams.moveEndTime ?
                vector<string> { to_string(sc.boundParams.moveEndTime) } :
                vector<string> {};
        }
    );
    boundParser.addStringArgsWithAliases(
        "TRANSFERSHAREAXIS", { "TRANSFERSHAREAXIS:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            if(lineVector.size() != 3) {
                cout <<
                     "There was an error parsing input file at Chemistry parameters. Exiting."
                     << endl;
                exit(EXIT_FAILURE);
            }

            else{
                cout<<"TRANSFERSHARE AXIS "<<lineVector[2]<<endl;
                if(lineVector[2]=="X")
                    sc.boundParams.transfershareaxis=0;
                else if(lineVector[2]=="Y")
                    sc.boundParams.transfershareaxis=1;
                else if(lineVector[2]=="Z")
                    sc.boundParams.transfershareaxis=2;
                else if(lineVector[2]=="RADIAL") {
                    sc.boundParams.transfershareaxis = 3;
                    cout<<"RADIAL transfer not implemented. Change paramters. Exiting"
                            "."<<endl;
                    exit(EXIT_FAILURE);
                }
                else{
                    cout << "There was an error parsing input file at Chemistry parameters. Exiting."
                         << endl;
                    exit(EXIT_FAILURE);
                }
            }
        },
        [] (const SimulConfig& sc) {
            vector<vector<string>> res;
            switch(sc.boundParams.transfershareaxis) {
                // TODO: figure out the 1st element
            case 0:
                res.push_back({ "???", "X" });
                break;
            case 1:
                res.push_back({ "???", "Y" });
                break;
            case 2:
                res.push_back({ "???", "Z" });
                break;
            case 3:
                res.push_back({ "???", "RADIAL" });
                break;
            }
            return res;
        }
    );
    boundParser.addEmptyLine();

    boundParser.addComment("# Set filament creation bounds.");
    boundParser.addComment("# Changes are not recommended.");
    boundParser.addEmptyLine();

    boundParser.addStringArgsWithAliases(
        "FILCREATIONBOUNDS", { "FILCREATIONBOUNDS:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            sc.boundParams.fraccompartmentspan.clear();
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
                sc.boundParams.fraccompartmentspan = tempbounds;
            }
        },
        [] (const SimulConfig& sc) {
            vector<vector<string>> res;
            if(const auto& frac = sc.boundParams.fraccompartmentspan; !frac.empty()) {
                res.push_back({
                    to_string(frac[0][0]), to_string(frac[0][1]), to_string(frac[0][2]),
                    to_string(frac[1][0]), to_string(frac[1][1]), to_string(frac[1][2])
                });
            }
            return res;
        }
    );
    boundParser.addEmptyLine();
}

void SystemParser::boundPostProcessing(SimulConfig& sc) const {

    if(const auto& bm = sc.boundParams.boundaryType.boundaryMove; !bm.empty()) {
        vector<int> leftfrontbottom = {0,0,0};
        vector<int> rightbacktop = {0,0,0};

        for(const auto& eachBM : bm) {
            if(eachBM == "LEFT")
                leftfrontbottom[0] = 1;
            else if(eachBM == "BOTTOM")
                leftfrontbottom[1] = 1;
            else if(eachBM == "FRONT")
                leftfrontbottom[2] = 1;
            else if(eachBM == "RIGHT")
                rightbacktop[0] = 1;
            else if(eachBM == "TOP")
                rightbacktop[1] = 1;
            else if(eachBM == "BACK")
                rightbacktop[2] = 1;
        }

        for(int i = 0; i < 3; i++){
            int addthemup = leftfrontbottom[i] + rightbacktop[i];
            if(addthemup > 0)
                sc.boundParams.transfershareaxis = i;
            if(addthemup == 2)
                sc.boundParams.planestomove = 2;
            else if(leftfrontbottom[i] == 1)
                sc.boundParams.planestomove = 1;
            else if(rightbacktop[i] == 1)
                sc.boundParams.planestomove = 0;
        }
    }

}

void SystemParser::initDyRateParser() {
    using namespace std;

    dyRateParser.addComment("##################################################");
    dyRateParser.addComment("### Dynamic rate parameters");
    dyRateParser.addComment("##################################################");
    dyRateParser.addEmptyLine();

    dyRateParser.addComment("# Brownian ratchet model");
    dyRateParser.addComment("# Footer et al, PNAS, 2007");
    dyRateParser.addEmptyLine();
    dyRateParser.addStringArgsWithAliases(
        "DFPOLYMERIZATIONTYPE", { "DFPOLYMERIZATIONTYPE:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            if (lineVector.size() >= 2) {
                for(int i = 1; i < lineVector.size(); i++)
                    sc.dyRateParams.dynamicRateType.dFPolymerizationType.push_back(lineVector[i]);
            }
        },
        [] (const SimulConfig& sc) {
            return sc.dyRateParams.dynamicRateType.dFPolymerizationType;
        }
    );
    dyRateParser.addStringArgsWithAliases(
        "DFPOLYMERIZATIONLEN", { "DFPOLYMERIZATIONLEN:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            if (lineVector.size() >= 2) {
                for(int i = 1; i < lineVector.size(); i++)
                    sc.dyRateParams.dFilPolymerizationCharLength.push_back(
                            atof((lineVector[i].c_str())));
            }
        },
        [] (const SimulConfig& sc) {
            vector<string> res;
            for(const auto& x : sc.dyRateParams.dFilPolymerizationCharLength) {
                res.push_back(to_string(x));
            }
            return res;
        }
    );
    dyRateParser.addEmptyLine();

    dyRateParser.addComment("# Motor catch-bond model");
    dyRateParser.addComment("# Erdmann et al, JCP, 2013");
    dyRateParser.addEmptyLine();
    dyRateParser.addStringArgsWithAliases(
        "DMUNBINDINGTYPE", { "DMUNBINDINGTYPE:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            if (lineVector.size() >= 2) {
                for(int i = 1; i < lineVector.size(); i++)
                    sc.dyRateParams.dynamicRateType.dMUnbindingType.push_back(lineVector[i]);
            }
        },
        [] (const SimulConfig& sc) {
            return sc.dyRateParams.dynamicRateType.dMUnbindingType;
        }
    );
    dyRateParser.addStringArgsWithAliases(
        "DMUNBINDINGFORCE", { "DMUNBINDINGFORCE:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            if (lineVector.size() >= 2) {
                for(int i = 1; i < lineVector.size(); i++)
                    sc.dyRateParams.dMotorUnbindingCharForce.push_back(
                            atof((lineVector[i].c_str())));
            }
        },
        [] (const SimulConfig& sc) {
            vector<string> res;
            for(const auto& x : sc.dyRateParams.dMotorUnbindingCharForce) {
                res.push_back(to_string(x));
            }
            return res;
        }
    );
    dyRateParser.addEmptyLine();

    dyRateParser.addComment("# Motor walking");
    dyRateParser.addComment("# Komianos & Papoian, PRX, 2018");
    dyRateParser.addComment("# Tunable parameters based on different studies");
    dyRateParser.addComment("# recommended values 24pN - 100 pN");
    dyRateParser.addEmptyLine();
    dyRateParser.addStringArgsWithAliases(
        "DMWALKINGTYPE", { "DMWALKINGTYPE:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            if (lineVector.size() >= 2) {
                for(int i = 1; i < lineVector.size(); i++)
                    sc.dyRateParams.dynamicRateType.dMWalkingType.push_back(lineVector[i]);
            }
        },
        [] (const SimulConfig& sc) {
            return sc.dyRateParams.dynamicRateType.dMWalkingType;
        }
    );
    dyRateParser.addStringArgsWithAliases(
        "DMWALKINGFORCE", { "DMWALKINGFORCE:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            if (lineVector.size() >= 2) {
                for(int i = 1; i < lineVector.size(); i++)
                    sc.dyRateParams.dMotorWalkingCharForce.push_back(
                            atof((lineVector[i].c_str())));
            }
        },
        [] (const SimulConfig& sc) {
            vector<string> res;
            for(const auto& x : sc.dyRateParams.dMotorWalkingCharForce) {
                res.push_back(to_string(x));
            }
            return res;
        }
    );
    dyRateParser.addEmptyLine();

    dyRateParser.addComment("# Linker slip-bond model");
    dyRateParser.addComment("# Ferrer et al, PNAS, 2008");
    dyRateParser.addEmptyLine();
    dyRateParser.addStringArgsWithAliases(
        "DLUNBINDINGTYPE", { "DLUNBINDINGTYPE:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            if (lineVector.size() >= 2) {
                for(int i = 1; i < lineVector.size(); i++)
                    sc.dyRateParams.dynamicRateType.dLUnbindingType.push_back(lineVector[i]);
            }
        },
        [] (const SimulConfig& sc) {
            return sc.dyRateParams.dynamicRateType.dLUnbindingType;
        }
    );
    dyRateParser.addStringArgsWithAliases(
        "DLUNBINDINGLEN", { "DLUNBINDINGLEN:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            if (lineVector.size() >= 2) {
                for(int i = 1; i < lineVector.size(); i++)
                    sc.dyRateParams.dLinkerUnbindingCharLength.push_back(
                            atof((lineVector[i].c_str())));
            }
        },
        [] (const SimulConfig& sc) {
            vector<string> res;
            for(const auto& x : sc.dyRateParams.dLinkerUnbindingCharLength) {
                res.push_back(to_string(x));
            }
            return res;
        }
    );
    dyRateParser.addStringArgsWithAliases(
        "DLUNBINDINGAMP", { "DLUNBINDINGAMP:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            if (lineVector.size() >= 2) {
                for(int i = 1; i < lineVector.size(); i++)
                    sc.dyRateParams.dLinkerUnbindingAmplitude.push_back(
                            atof((lineVector[i].c_str())));
            }
        },
        [] (const SimulConfig& sc) {
            vector<string> res;
            for(const auto& x : sc.dyRateParams.dLinkerUnbindingAmplitude) {
                res.push_back(to_string(x));
            }
            return res;
        }
    );
    dyRateParser.addEmptyLine();

    dyRateParser.addComment("# Brancher slip-bond model");
    dyRateParser.addComment("# Source unclear.");
    dyRateParser.addEmptyLine();
    dyRateParser.addStringArgsWithAliases(
        "DBUNBINDINGTYPE", { "DBUNBINDINGTYPE:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            if (lineVector.size() >= 2) {
                for(int i = 1; i < lineVector.size(); i++)
                    sc.dyRateParams.dynamicRateType.dBUnbindingType.push_back(lineVector[i]);
            }
        },
        [] (const SimulConfig& sc) {
            return sc.dyRateParams.dynamicRateType.dBUnbindingType;
        }
    );
    dyRateParser.addStringArgsWithAliases(
        "DBUNBINDINGLEN", { "DBUNBINDINGLEN:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            if (lineVector.size() >= 2) {
                for(int i = 1; i < lineVector.size(); i++)
                    sc.dyRateParams.dBranchUnbindingCharLength.push_back(
                            atof((lineVector[i].c_str())));
            }
        },
        [] (const SimulConfig& sc) {
            vector<string> res;
            for(const auto& x : sc.dyRateParams.dBranchUnbindingCharLength) {
                res.push_back(to_string(x));
            }
            return res;
        }
    );
    dyRateParser.addStringArgsWithAliases(
        "DBUNBINDINGF", { "DBUNBINDINGF:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            if (lineVector.size() >= 2) {
                for(int i = 1; i < lineVector.size(); i++)
                    sc.dyRateParams.dBranchUnbindingCharForce.push_back(
                            atof((lineVector[i].c_str())));
            }
        },
        [] (const SimulConfig& sc) {
            vector<string> res;
            for(const auto& x : sc.dyRateParams.dBranchUnbindingCharForce) {
                res.push_back(to_string(x));
            }
            return res;
        }
    );
    dyRateParser.addEmptyLine();

    /// Manual Rate Changer
    // It takes 5 inputs as start_time, plusend_poly, plusend_depoly, minusend_poly, minusend_depoly
    // Currently it applies type 0 to all filament types
    dyRateParser.addComment("### Manual rate changer.");
    dyRateParser.addComment("# Do not change unless you know what you are doing.");
    dyRateParser.addEmptyLine();
    dyRateParser.addStringArgsWithAliases(
        "MANUALSTARTTIME", { "MANUALSTARTTIME:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            if (lineVector.size() >= 2) {
                sc.dyRateParams.manualCharStartTime = atof((lineVector[1].c_str()));
            }
        },
        [] (const SimulConfig& sc) {
            return vector<string> { to_string(sc.dyRateParams.manualCharStartTime) };
        }
    );
    dyRateParser.addStringArgsWithAliases(
        "MANUALPLUSPOLYRATIO", { "MANUALPLUSPOLYRATIO:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            if (lineVector.size() >= 2) {
                sc.dyRateParams.manualPlusPolyRate = atof((lineVector[1].c_str()));
            }
        },
        [] (const SimulConfig& sc) {
            return vector<string> { to_string(sc.dyRateParams.manualPlusPolyRate) };
        }
    );
    dyRateParser.addStringArgsWithAliases(
        "MANUALPLUSDEPOLYRATIO", { "MANUALPLUSDEPOLYRATIO:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            if (lineVector.size() >= 2) {
                sc.dyRateParams.manualPlusDepolyRate = atof((lineVector[1].c_str()));
            }
        },
        [] (const SimulConfig& sc) {
            return vector<string> { to_string(sc.dyRateParams.manualPlusDepolyRate) };
        }
    );
    dyRateParser.addStringArgsWithAliases(
        "MANUALMINUSPOLYRATIO", { "MANUALMINUSPOLYRATIO:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            if (lineVector.size() >= 2) {
                sc.dyRateParams.manualMinusPolyRate = atof((lineVector[1].c_str()));
            }
        },
        [] (const SimulConfig& sc) {
            return vector<string> { to_string(sc.dyRateParams.manualMinusPolyRate) };
        }
    );
    dyRateParser.addStringArgsWithAliases(
        "MANUALMINUSDEPOLYRATIO", { "MANUALMINUSDEPOLYRATIO:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            if (lineVector.size() >= 2) {
                sc.dyRateParams.manualMinusDepolyRate = atof((lineVector[1].c_str()));
            }
        },
        [] (const SimulConfig& sc) {
            return vector<string> { to_string(sc.dyRateParams.manualMinusDepolyRate) };
        }
    );
    dyRateParser.addEmptyLine();

}

void SystemParser::initGeoParser() {
    using namespace std;

    geoParser.addComment("##################################################");
    geoParser.addComment("### Geometric parameters");
    geoParser.addComment("##################################################");
    geoParser.addEmptyLine();

    geoParser.addComment("### Set network sizes and shape");
    geoParser.addComment("# Set the number of compartments in x, y and z directions");
    geoParser.addComment("# Network size = compartment size (500nm by default) * (NX, NY, NZ)");
    geoParser.addEmptyLine();

    geoParser.addStringArgsWithAliases(
        "NX", { "NX:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            if(lineVector.size() > 2) {
                cout << "There was an error parsing input file at grid dimensions. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
            else if(lineVector.size() == 2)
                sc.geoParams.NX = stoi(lineVector[1]);
            else {}
        },
        [] (const SimulConfig& sc) {
            vector< string > res;
            if(sc.geoParams.NX) res.push_back(to_string(sc.geoParams.NX));
            return res;
        }
    );
    geoParser.addStringArgsWithAliases(
        "NY", { "NY:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            if(lineVector.size() > 2) {
                cout << "There was an error parsing input file at grid dimensions. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
            else if(lineVector.size() == 2)
                sc.geoParams.NY = stoi(lineVector[1]);
            else {}
        },
        [] (const SimulConfig& sc) {
            vector< string > res;
            if(sc.geoParams.NY) res.push_back(to_string(sc.geoParams.NY));
            return res;
        }
    );
    geoParser.addStringArgsWithAliases(
        "NZ", { "NZ:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            if(lineVector.size() > 2) {
                cout << "There was an error parsing input file at grid dimensions. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
            else if(lineVector.size() == 2)
                sc.geoParams.NZ = stoi(lineVector[1]);
            else {}
        },
        [] (const SimulConfig& sc) {
            vector< string > res;
            if(sc.geoParams.NZ) res.push_back(to_string(sc.geoParams.NZ));
            return res;
        }
    );
    geoParser.addEmptyLine();

    geoParser.addComment("### The compartment size");
    geoParser.addComment("# Based on Kuramoto length, see Popov et al., PLoS Comp Biol, 2016 ");
    geoParser.addComment("# Some chemical reaction rates are scaled based on compartment size");
    geoParser.addEmptyLine();

    geoParser.addStringArgsWithAliases(
        "COMPARTMENTSIZEX", { "COMPARTMENTSIZEX:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            if(lineVector.size() > 2) {
                cout << "There was an error parsing input file at compartment size. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
            else if(lineVector.size() == 2)
                sc.geoParams.compartmentSizeX = stod(lineVector[1]);
            else {}
        },
        [] (const SimulConfig& sc) {
            vector< string > res;
            if(sc.geoParams.compartmentSizeX) res.push_back(to_string(sc.geoParams.compartmentSizeX));
            return res;
        }
    );
    geoParser.addStringArgsWithAliases(
        "COMPARTMENTSIZEY", { "COMPARTMENTSIZEY:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            if(lineVector.size() > 2) {
                cout << "There was an error parsing input file at compartment size. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
            else if(lineVector.size() == 2)
                sc.geoParams.compartmentSizeY = stod(lineVector[1]);
            else {}
        },
        [] (const SimulConfig& sc) {
            vector< string > res;
            if(sc.geoParams.compartmentSizeY) res.push_back(to_string(sc.geoParams.compartmentSizeY));
            return res;
        }
    );
    geoParser.addStringArgsWithAliases(
        "COMPARTMENTSIZEZ", { "COMPARTMENTSIZEZ:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            if(lineVector.size() > 2) {
                cout << "There was an error parsing input file at compartment size. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
            else if(lineVector.size() == 2)
                sc.geoParams.compartmentSizeZ = stod(lineVector[1]);
            else {}
        },
        [] (const SimulConfig& sc) {
            vector< string > res;
            if(sc.geoParams.compartmentSizeZ) res.push_back(to_string(sc.geoParams.compartmentSizeZ));
            return res;
        }
    );
    geoParser.addEmptyLine();

    geoParser.addComment("### Cylinder setup");
    geoParser.addComment("### Changes not recommended");
    geoParser.addEmptyLine();

    geoParser.addStringArgsWithAliases(
        "MONOMERSIZE", { "MONOMERSIZE:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            sc.geoParams.monomerSize.clear();
            if (lineVector.size() >= 2) {
                for(int i = 1; i < lineVector.size(); i++)
                    sc.geoParams.monomerSize.push_back(atof((lineVector[i].c_str())));
            }
            else {}
        },
        [] (const SimulConfig& sc) {
            vector< string > res;
            for(const auto& s : sc.geoParams.monomerSize)
                res.push_back(to_string(s));
            return res;
        }
    );
    geoParser.addStringArgsWithAliases(
        "CYLINDERSIZE", { "CYLINDERSIZE:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            sc.geoParams.cylinderSize.clear();
            if (lineVector.size() >= 2) {
                for(int i = 1; i < lineVector.size(); i++)
                    sc.geoParams.cylinderSize.push_back(atof((lineVector[i].c_str())));
            }
            else {}
        },
        [] (const SimulConfig& sc) {
            vector< string > res;
            for(const auto& s : sc.geoParams.cylinderSize)
                res.push_back(to_string(s));
            return res;
        }
    );
    geoParser.addEmptyLine();

    geoParser.addComment("### Simulation dimension");
    geoParser.addComment("### DO NOT CHANGE");
    geoParser.addEmptyLine();

    geoParser.addStringArgsWithAliases(
        "NDIM", { "NDIM:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            if(lineVector.size() !=  2) {
                cout << "Number of dimensions needs to be specified. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
            else {
                sc.geoParams.nDim = short(atoi(lineVector[1].c_str()));
            }
        },
        [] (const SimulConfig& sc) {
            return vector< string > { to_string(sc.geoParams.nDim) };
        }
    );
    geoParser.addEmptyLine();
}

void SystemParser::geoPostProcessing(SimulConfig& sc) const {
    using namespace std;

    if(sc.geoParams.cylinderSize.size() != sc.geoParams.monomerSize.size()) {

        cout << "Must specify an equivalent number of cylinder and monomer sizes. Exiting."
             << endl;
        exit(EXIT_FAILURE);
    }

    for(int i = 0; i < sc.geoParams.cylinderSize.size(); i++) {

#ifdef CHEMISTRY
        if(sc.geoParams.cylinderSize[i] / sc.geoParams.monomerSize[i] < sc.geoParams.minCylinderNumMon) {
            cout <<
                 "With chemistry, cylinder size specified is too short. Exiting."
                 << endl;
            exit(EXIT_FAILURE);
        }
#endif
        sc.geoParams.cylinderNumMon.push_back(int(sc.geoParams.cylinderSize[i] / sc.geoParams.monomerSize[i]));

        sc.geoParams.minCylinderSize.push_back(
                sc.geoParams.minCylinderNumMon * sc.geoParams.monomerSize[i]);

    }

    //find max compartment side
    sc.geoParams.largestCompartmentSide = max({
        sc.geoParams.compartmentSizeX,
        sc.geoParams.compartmentSizeY,
        sc.geoParams.compartmentSizeZ
    });
    //find max Cylinder size
    sc.geoParams.largestCylinderSize = 0;
    for(auto l : sc.geoParams.cylinderSize)
        sc.geoParams.largestCylinderSize = max(sc.geoParams.largestCylinderSize, l);

}

void SystemParser::initInitParser() {
    initParser.addComment("##################################################");
    initParser.addComment("### Initial network setup");
    initParser.addComment("##################################################");
    initParser.addEmptyLine();

    initParser.addComment("# Initial bubble setup from external input");
    initParser.addStringArgsWithAliases(
        "BUBBLEFILE", { "BUBBLEFILE:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            if(lineVector.size() > 2) {
                cout << "Error reading bubble input file. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2)
                sc.bubbleSetup.inputFile = lineVector[1];
            else {}
        },
        [] (const SimulConfig& sc) {
            vector<string> res;
            if(!sc.bubbleSetup.inputFile.empty()) {
                res.push_back( sc.bubbleSetup.inputFile.string() );
            }
            return res;
        }
    );
    initParser.addEmptyLine();

    initParser.addComment("# Generate initial bubbles");
    initParser.addStringArgsWithAliases(
        "NUMBUBBLES", { "NUMBUBBLES:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            if(lineVector.size() > 2) {
                cout << "Error reading number of bubbles. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2)
                sc.bubbleSetup.numBubbles = atoi(lineVector[1].c_str());
            else {}
        },
        [] (const SimulConfig& sc) {
            vector<string> res;
            if(sc.bubbleSetup.numBubbles) {
                res.push_back( to_string(sc.bubbleSetup.numBubbles) );
            }
            return res;
        }
    );
    initParser.addStringArgsWithAliases(
        "BUBBLETYPE", { "BUBBLETYPE:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            if(lineVector.size() > 2) {
                cout << "Error reading bubble type. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2)
                sc.bubbleSetup.bubbleType = atoi(lineVector[1].c_str());
            else {}
        },
        [] (const SimulConfig& sc) {
            vector<string> res;
            if(sc.bubbleSetup.bubbleType) {
                res.push_back( to_string(sc.bubbleSetup.bubbleType) );
            }
            return res;
        }
    );
    initParser.addEmptyLine();

    initParser.addComment("# Initial filament setup from external input");
    initParser.addStringArgsWithAliases(
        "FILAMENTFILE", { "FILAMENTFILE:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            if(lineVector.size() > 2) {
                cout << "Error reading filament input file. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2)
                sc.filamentSetup.inputFile = lineVector[1];
            else {}
        },
        [] (const SimulConfig& sc) {
            vector<string> res;
            if(!sc.filamentSetup.inputFile.empty()) {
                res.push_back( sc.filamentSetup.inputFile.string() );
            }
            return res;
        }
    );
    initParser.addEmptyLine();

    initParser.addComment("# Generate initial filaments");
    initParser.addStringArgsWithAliases(
        "NUMFILAMENTS", { "NUMFILAMENTS:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            if(lineVector.size() > 2) {
                cout << "Error reading number of filaments. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2)
                sc.filamentSetup.numFilaments = atoi(lineVector[1].c_str());
            else {}
        },
        [] (const SimulConfig& sc) {
            vector<string> res;
            if(sc.filamentSetup.numFilaments) {
                res.push_back( to_string(sc.filamentSetup.numFilaments) );
            }
            return res;
        }
    );
    initParser.addStringArgsWithAliases(
        "FILAMENTLENGTH", { "FILAMENTLENGTH:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            if(lineVector.size() > 2) {
                cout << "Error reading filament length. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2)
                sc.filamentSetup.filamentLength = atoi(lineVector[1].c_str());
            else {}
        },
        [] (const SimulConfig& sc) {
            vector<string> res;
            if(sc.filamentSetup.filamentLength) {
                res.push_back( to_string(sc.filamentSetup.filamentLength) );
            }
            return res;
        }
    );
    initParser.addStringArgsWithAliases(
        "FILAMENTTYPE", { "FILAMENTTYPE:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            if(lineVector.size() > 2) {
                cout << "Error reading filament type. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2)
                sc.filamentSetup.filamentType = atoi(lineVector[1].c_str());
            else {}
        },
        [] (const SimulConfig& sc) {
            vector<string> res;
            if(sc.filamentSetup.filamentType) {
                res.push_back( to_string(sc.filamentSetup.filamentType) );
            }
            return res;
        }
    );
    initParser.addStringArgsWithAliases(
        "PROJECTIONTYPE", { "PROJECTIONTYPE:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            if(lineVector.size() > 2) {
                cout << "Error reading filament projection type. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2)
                sc.filamentSetup.projectionType = lineVector[1];
            else {}
        },
        [] (const SimulConfig& sc) {
            vector<string> res;
            if(!sc.filamentSetup.projectionType.empty()) {
                res.push_back( sc.filamentSetup.projectionType );
            }
            return res;
        }
    );
    initParser.addEmptyLine();

    initParser.addComment("# Restart settings");
    initParser.addStringArgsWithAliases(
        "RESTARTPHASE", { "RESTARTPHASE:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            if(lineVector.size() > 2) {
                cout << "Error reading restart params. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2){
                if(lineVector[1].find("USECHEMCOPYNUM"))
                sc.filamentSetup.USECHEMCOPYNUM = true;
            }
        },
        [] (const SimulConfig& sc) {
            vector<string> res;
            if(sc.filamentSetup.USECHEMCOPYNUM) {
                res.push_back("USECHEMCOPYNUM");
            }
            return res;
        }
    );
    initParser.addStringArgsWithAliases(
        "PINRESTARTFILE", { "PINRESTARTFILE:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            if(lineVector.size() > 2) {
                cout << "Error reading filament projection type. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2)
                sc.filamentSetup.pinRestartFile = lineVector[1];
            else {}
        },
        [] (const SimulConfig& sc) {
            vector<string> res;
            if(!sc.filamentSetup.pinRestartFile.empty()) {
                res.push_back(sc.filamentSetup.pinRestartFile);
            }
            return res;
        }
    );
    initParser.addEmptyLine();

    initParser.addComment("# Special setup types");
    initParser.addStringArgsWithAliases(
        "SPECIALSETUP", { "SPECIALSETUP:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            if(lineVector.size() > 3) {
                cout <<
                     "There was an error parsing input file at special setup types. Exiting."
                     << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector[1] == "MTOC") sc.specialParams.specialSetupType.mtoc = true;
            else if (lineVector[1] == "AFM") sc.specialParams.specialSetupType.afm = true;
        },
        [] (const SimulConfig& sc) {
            vector<vector<string>> res;
            if(sc.specialParams.specialSetupType.mtoc) res.push_back({ "MTOC" });
            if(sc.specialParams.specialSetupType.afm) res.push_back({ "AFM" });
            return res;
        }
    );
    initParser.addStringArgsWithAliases(
        "MTOCFILAMENTTYPE", { "MTOCFILAMENTTYPE:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            if(lineVector.size() != 2) {
                cout <<
                     "A filament type to connect to an MTOC must be specified. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                sc.specialParams.specialSetupType.mtocFilamentType = atoi(lineVector[1].c_str());
            }
        },
        [] (const SimulConfig& sc) {
            return vector<string> { to_string(sc.specialParams.specialSetupType.mtocFilamentType) };
        }
    );
    initParser.addStringArgsWithAliases(
        "MTOCNUMFILAMENTS", { "MTOCNUMFILAMENTS:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            if(lineVector.size() != 2) {
                cout <<
                     "A number of filaments to connect to an MTOC must be specified. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                sc.specialParams.specialSetupType.mtocNumFilaments = atoi(lineVector[1].c_str());
            }
        },
        [] (const SimulConfig& sc) {
            return vector<string> { to_string(sc.specialParams.specialSetupType.mtocNumFilaments) };
        }
    );
    initParser.addStringArgsWithAliases(
        "MTOCFILAMENTLENGTH", { "MTOCFILAMENTLENGTH:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            if(lineVector.size() != 2) {
                cout <<
                     "A filament length for filaments connected to an MTOC must be specified. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                sc.specialParams.specialSetupType.mtocFilamentLength = atoi(lineVector[1].c_str());
            }
        },
        [] (const SimulConfig& sc) {
            return vector<string> { to_string(sc.specialParams.specialSetupType.mtocFilamentLength) };
        }
    );
    initParser.addStringArgsWithAliases(
        "MTOCBUBBLETYPE", { "MTOCBUBBLETYPE:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            if(lineVector.size() != 2) {
                cout <<
                     "A bubble type to connect to an MTOC must be specified. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                sc.specialParams.specialSetupType.mtocBubbleType = atoi(lineVector[1].c_str());
            }
        },
        [] (const SimulConfig& sc) {
            return vector<string> { to_string(sc.specialParams.specialSetupType.mtocBubbleType) };
        }
    );
    initParser.addStringArgsWithAliases(
        "AFMFILAMENTTYPE", { "AFMFILAMENTTYPE:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            if(lineVector.size() != 2) {
                cout <<
                    "A filament type to connect to an AFM must be specified. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                sc.specialParams.specialSetupType.afmFilamentType = atoi(lineVector[1].c_str());
            }
        },
        [] (const SimulConfig& sc) {
            return vector<string> { to_string(sc.specialParams.specialSetupType.afmFilamentType) };
        }
    );
    initParser.addStringArgsWithAliases(
        "AFMNUMFILAMENTS", { "AFMNUMFILAMENTS:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            if(lineVector.size() != 2) {
                cout <<
                "A number of filaments to connect to an AFM must be specified. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                sc.specialParams.specialSetupType.afmNumFilaments = atoi(lineVector[1].c_str());
            }
        },
        [] (const SimulConfig& sc) {
            return vector<string> { to_string(sc.specialParams.specialSetupType.afmNumFilaments) };
        }
    );
    initParser.addStringArgsWithAliases(
        "AFMFILAMENTLENGTH", { "AFMFILAMENTLENGTH:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            if(lineVector.size() != 2) {
                cout <<
                "A filament length for filaments connected to an AFM must be specified. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                sc.specialParams.specialSetupType.afmFilamentLength = atoi(lineVector[1].c_str());
            }
        },
        [] (const SimulConfig& sc) {
            return vector<string> { to_string(sc.specialParams.specialSetupType.afmFilamentLength) };
        }
    );
    initParser.addStringArgsWithAliases(
        "AFMBUBBLETYPE", { "AFMBUBBLETYPE:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            if(lineVector.size() != 2) {
                cout <<
                "A bubble type to connect to an MTOC must be specified. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 2) {
                sc.specialParams.specialSetupType.afmBubbleType = atoi(lineVector[1].c_str());
            }
        },
        [] (const SimulConfig& sc) {
            return vector<string> { to_string(sc.specialParams.specialSetupType.afmBubbleType) };
        }
    );
    initParser.addStringArgsWithAliases(
        "MTOCXYZCOORD", { "MTOCXYZCOORD:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            sc.specialParams.specialSetupType.mtocInputCoordXYZ.clear();
            if(lineVector.size() != 4) {
                cout <<
                "Coordinates of MTOC must be specified. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 4) {
                sc.specialParams.specialSetupType.mtocInputCoordXYZ.push_back(stof(lineVector[1].c_str()));
                sc.specialParams.specialSetupType.mtocInputCoordXYZ.push_back(stof(lineVector[2].c_str()));
                sc.specialParams.specialSetupType.mtocInputCoordXYZ.push_back(stof(lineVector[3].c_str()));
            }
        },
        [] (const SimulConfig& sc) {
            vector<vector<string>> res;
            if(const auto& xyz = sc.specialParams.specialSetupType.mtocInputCoordXYZ; !xyz.empty()) {
                res.push_back({
                    to_string(xyz[0]),
                    to_string(xyz[1]),
                    to_string(xyz[2])
                });
            }
            return res;
        }
    );
    initParser.addEmptyLine();

    initParser.addComment("# Special setup parameters");
    initParser.addStringArgsWithAliases(
        "MTOCFILAMENTCOORD", { "MTOCFILAMENTCOORD:" },
        [] (SimulConfig& sc, const vector<string>& lineVector) {
            if(lineVector.size() != 5) {
                cout << "4 coordinates of MTOC filaments must be specified. Exiting." << endl;
                exit(EXIT_FAILURE);
            }
            else if (lineVector.size() == 5) {
                sc.specialParams.mtocTheta1 = stof(lineVector[1].c_str());
                sc.specialParams.mtocTheta2 = stof(lineVector[2].c_str());
                sc.specialParams.mtocPhi1 = stof(lineVector[3].c_str());
                sc.specialParams.mtocPhi2 = stof(lineVector[4].c_str());
            }
        },
        [] (const SimulConfig& sc) {
            return vector<string> {
                to_string(sc.specialParams.mtocTheta1),
                to_string(sc.specialParams.mtocTheta2),
                to_string(sc.specialParams.mtocPhi1),
                to_string(sc.specialParams.mtocPhi2)
            };
        }
    );

}

} // namespace medyan

FilamentData FilamentParser::readFilaments(std::istream& is) {
    is.clear();
    is.seekg(0);
     vector<tuple<short, vector<floatingpoint>, vector<floatingpoint>>> filamentVector;
     vector<vector<vector<floatingpoint>>> linkerVector;
     vector<vector<vector<floatingpoint>>> motorVector;
     vector<vector<floatingpoint>> staticVector;
     vector<tuple<string, short, vector<vector<floatingpoint>>>> boundVector;
     vector<tuple<string, short, vector<floatingpoint>>> branchVector;
     string line;
    
    while(getline(is, line)) {

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

vector<tuple<short, vector<floatingpoint>>> BubbleParser::readBubbles(std::istream& is) {
    
    is.clear();
    is.seekg(0);
    
    vector<tuple<short, vector<floatingpoint>>> returnVector;
    string line;

    while(getline(is, line)) {

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

ChemistryData ChemistryParser::readChemistryInput(std::istream& is, const ChemParams& chemParams) {

    is.clear();
    is.seekg(0);

    ///To keep track of duplicate names
    vector<string> allSpeciesNames;

    ChemistryData chem; string line;

    while(getline(is, line)) {

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
            if(chemParams.dissTracking){
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
            if(chemParams.dissTracking){
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
            if(chemParams.dissTracking){
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
            if(chemParams.dissTracking){
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
            if(chemParams.dissTracking){
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
            if(chemParams.dissTracking){
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
            if(chemParams.dissTracking){
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
