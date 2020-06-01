#ifndef MEDYAN_MedyanConfig_hpp
#define MEDYAN_MedyanConfig_hpp

/*
This file includes all one needs for interacting with MEDYAN user inputs.

Features:
  - Simulation config storage (SysParams.h)
  - Simulation config validation (SysParams.h)
  - Reading config from input files (Parser.h)
  - Reading config from command line guide
  - Generating config input files
  - GUI config display (future)
  - GUI config interactive modification (very future)

Documentation:

  Data structures
  ---------------
    medyan::SimulConfig

      All the simulation configuration is stored in an object of this type.

    medyan::SimulConfigHelper

      Manages the parsers and does configuration input/output.

*/

#include <filesystem>
#include <iostream>
#include <stdexcept>

#include "Parser.h"
#include "SysParams.h"
#include "Util/Io/Log.hpp"

namespace medyan {


//---------------------------------
// Read from input
//---------------------------------

// Auxiliary: read the whole file into string
inline std::string readFileToString(std::filesystem::path file) {
    using namespace std;

    ifstream ifs(file);
    if(!ifs.is_open()) {
        LOG(ERROR) << "There was an error reading file " << file;
        throw std::runtime_error("Cannot open input file.");
    }
    LOG(INFO) << "Loading file " << file;

    stringstream ss;
    ss << ifs.rdbuf(); // Read the file
    return ss.str();
}


// Read bubble input
inline void readBubbleConfig(SimulConfig& sc, std::istream& is) {
    sc.bubbleData = BubbleParser::readBubbles(is);
}

// Read filament input
inline void readFilamentConfig(SimulConfig& sc, std::istream& is) {
    sc.filamentData = FilamentParser::readFilaments(is);
}

struct SimulConfigHelper {

    // Parsers
    SystemParser systemParser;
    ChemistryParser chemistryParser;

    // Read simulation configuration from file
    // Not const function because chemistry parsing is non-const
    SimulConfig getFromInput(
        std::filesystem::path systemInputFile,
        std::filesystem::path inputDirectory
    ) {
        using namespace std;
        using namespace std::filesystem;

        // Auxiliary struct handling file streams
        struct ReadFile {
            ifstream ifs;
            ReadFile(const path& p) : ifs(p) {
                if(!ifs.is_open()) {
                    LOG(ERROR) << "There was an error parsing file " << p;
                    throw std::runtime_error("Cannot open input file.");
                }
                LOG(INFO) << "Loading input file " << p;
            }
        };

        SimulConfig conf;
        conf.metaParams.systemInputFile = systemInputFile;
        conf.metaParams.inputDirectory  = inputDirectory;

        // Read system input
        {
            systemParser.parseInput(conf, readFileToString(systemInputFile));
        }

        // Read chemistry input
        if(conf.chemParams.chemistrySetup.inputFile.empty()) {
            LOG(ERROR) << "Need to specify a chemical input file. Exiting.";
            throw std::runtime_error("No chemistry input file specified.");
        }
        else {
            chemistryParser.parseInput(
                conf,
                readFileToString(inputDirectory / conf.chemParams.chemistrySetup.inputFile)
            );
        }

        // Read bubble input
        if(!conf.bubbleSetup.inputFile.empty()) {
            ReadFile file(inputDirectory / conf.bubbleSetup.inputFile);
            readBubbleConfig(conf, file.ifs);
        }

        // Read filament input
        if(!conf.filamentSetup.inputFile.empty()) {
            ReadFile file(inputDirectory / conf.filamentSetup.inputFile);
            readFilamentConfig(conf, file.ifs);
        }

        return conf;
    }

    // Output configuration to file
    void generateInput(const SimulConfig& conf) const {
        using namespace std;

        // Generate system input
        if(conf.metaParams.systemInputFile.empty()) {
            LOG(NOTE) << "The system input file is not specified in the config."
                " Using cout as the default output target.";
            systemParser.outputInput(cout, conf);
            LOG(NOTE) << "----- End of system input -----";
        }
        else {
            const auto& p = conf.metaParams.systemInputFile;
            if(filesystem::exists(p)) {
                LOG(WARNING) << "The file " << p << " already exists and will be overwritten.";
            }
            ofstream ofs(p);
            systemParser.outputInput(ofs, conf);
        }

        // Generate chemistry input
        if(conf.chemParams.chemistrySetup.inputFile.empty()) {
            LOG(NOTE) << "The chemistry input file is not specified in the config."
                " Using cout as the default output target.";
            chemistryParser.outputInput(cout, conf);
            LOG(NOTE) << "----- End of chemistry input -----";
        }
        else {
            const auto p = conf.metaParams.inputDirectory / conf.chemParams.chemistrySetup.inputFile;
            if(filesystem::exists(p)) {
                LOG(WARNING) << "The file " << p << " already exists and will be overwritten.";
            }
            ofstream ofs(p);
            chemistryParser.outputInput(ofs, conf);
        }

        // Generate bubble input
        if(!conf.bubbleSetup.inputFile.empty()) {
            LOG(WARNING) << "Additional bubble input is specified, but an input file is not generated.";
            LOG(INFO) << "It will be supported in future versions.";
        }

        // Generate filament input
        if(!conf.filamentSetup.inputFile.empty()) {
            LOG(WARNING) << "Additional filament input is specified, but an input file is not generate.";
            LOG(INFO) << "It will be supported in future versions.";
        }


    }

};


// Starts an interactive session of the command line configuration generation
inline void interactiveConfig() {
    using namespace std;

    // Auxiliary functions
    //-------------------------------------------------------------------------
    const auto getInt = [](int def) -> int {
        string line; getline(cin, line);
        if(line.empty()) return def;
        return stoi(line);
    };
    const auto getDouble = [](double def) -> double {
        string line; getline(cin, line);
        if(line.empty()) return def;
        return stod(line);
    };
    const auto getString = [](string def) -> string {
        string line; getline(cin, line);
        if(line.empty()) return def;
        return line;
    };
    const auto chooseStringAmong = [](vector<string> choices) -> string {
        if(choices.empty()) return "";
        // The first one is default
        auto s = getString(choices[0]);
        while(find(choices.begin(), choices.end(), s) == choices.end()) {
            cout << "Invalid choice: " << s << ". Please correct the input: ";
            s = getString(choices[0]);
        }
    };
    const auto chooseYesNo = [](bool def) -> bool {
        string line;
        getline(cin, line);
        if(line.empty()) return def;
        while(line != "y" && line != "n") {
            cout << "Please indicate yes or no using (y/n): ";
            getline(cin, line);
            if(line.empty()) return def;
        }
        return line == "y";
    };

    // States
    //-------------------------------------------------------------------------
    SimulConfigHelper helper;

    bool shouldAddMyosin = true;
    bool shouldAddLinker = true;
    bool shouldAddBrancher = true;

    filesystem::path fileDirectory;

    // Initialize a default configuration
    //-------------------------------------------------------------------------
    SimulConfig conf;

    conf.geoParams.nDim = 3;
    conf.geoParams.compartmentSizeX = 500.0;
    conf.geoParams.compartmentSizeY = 500.0;
    conf.geoParams.compartmentSizeZ = 500.0;
    conf.geoParams.monomerSize = {2.7};
    conf.geoParams.cylinderSize = {108.0};

    conf.chemParams.motorNumHeadsMin = {15};
    conf.chemParams.motorNumHeadsMax = {30};
    conf.chemParams.motorStepSize = {6.0};
    conf.chemParams.numBindingSites = {4};
    conf.chemParams.numFilaments = 1;
    conf.chemParams.chemistryAlgorithm.algorithm = "NRM";

    conf.mechParams.mechanicsFFType.BoundaryFFType = "REPULSIONEXP";
    conf.boundParams.BoundaryCutoff = 300.0;
    conf.boundParams.BoundaryK = 41.0;
    conf.boundParams.BScreenLength = 2.7;

    conf.mechParams.mechanicsFFType.FStretchingType = "HARMONIC";
    conf.mechParams.FStretchingK = {100.0};
    conf.mechParams.mechanicsFFType.FBendingType = "COSINE";
    conf.mechParams.FBendingK = {672.0};
    conf.mechParams.FBendingTheta = {0.0};

    conf.mechParams.mechanicsAlgorithm.ConjugateGradient = "POLAKRIBIERE";
    conf.mechParams.mechanicsAlgorithm.gradientTolerance = 5.0;
    conf.mechParams.mechanicsAlgorithm.maxDistance = 1.0;
    conf.mechParams.mechanicsAlgorithm.lambdaMax = 0.01;

    conf.mechParams.sameFilBindSkip = 100;

    conf.filamentSetup.projectionType = "STRAIGHT";

    conf.chemistryData.speciesFilament[0].push_back("AF");
    conf.chemistryData.speciesPlusEnd[0].push_back("PA");
    conf.chemistryData.speciesMinusEnd[0].push_back("MA");

    const auto addMyosin = [&] {
        conf.mechParams.mechanicsFFType.MStretchingType = "HARMONIC";
        conf.mechParams.MStretchingK = {2.5};

        conf.chemistryData.speciesMotor[0].push_back("MOA");
        conf.chemistryData.speciesBound[0].push_back("MEA");
        conf.chemistryData.M_BINDING_INDEX[0] = "MEA";

        conf.chemistryData.motorReactions.push_back({
            { "MEA:BOUND:1", "MEA:BOUND:2", "MD:DIFFUSING" },
            { "MOA:MOTOR:1", "MOA:MOTOR:2" },
            0.2, 1.7, 175.0, 225.0,
            0.0, "NA"
        });
        conf.chemistryData.motorWalkingReactions.push_back({
            { "MOA:MOTOR:N", "MEA:BOUND:N+1" },
            { "MOA:MOTOR:N+1", "MEA:BOUND:N" },
            0.2,
            0.0, "NA"
        });
    };
    const auto addLinker = [&] {
        conf.mechParams.mechanicsFFType.LStretchingType = "HARMONIC";
        conf.mechParams.LStretchingK = {8.0};

        conf.chemistryData.speciesLinker[0].push_back("LA");
        conf.chemistryData.speciesBound[0].push_back("LEA");
        conf.chemistryData.L_BINDING_INDEX[0] = "LEA";

        conf.chemistryData.linkerReactions.push_back({
            { "LEA:BOUND:1", "LEA:BOUND:2", "LD:DIFFUSING" },
            { "LA:LINKER:1", "LA:LINKER:2" },
            0.01, 0.3, 30.0, 40.0,
            0.0, "NA"
        });
    };
    const auto addBrancher = [&] {
        conf.mechParams.mechanicsFFType.BrStretchingType = "HARMONIC";
        conf.mechParams.BrStretchingK = {100.0};
        conf.mechParams.BrStretchingL = {6.0};
        conf.mechParams.mechanicsFFType.BrBendingType = "COSINE";
        conf.mechParams.BrBendingK = {10.0};
        conf.mechParams.BrBendingTheta = {1.22};
        conf.mechParams.mechanicsFFType.BrDihedralType = "COSINE";
        conf.mechParams.BrDihedralK = {10.0};
        conf.mechParams.mechanicsFFType.BrPositionType = "COSINE";
        conf.mechParams.BrPositionK = {20.0};

        conf.chemistryData.speciesBrancher[0].push_back("BA");
        conf.chemistryData.speciesBound[0].push_back("BEA");
        conf.chemistryData.B_BINDING_INDEX[0] = "BEA";

        conf.chemistryData.branchingReactions.push_back({
            { "BD:DIFFUSING", "A:DIFFUSING", "BEA:BOUND" },
            { "BA:BRANCHER", "PA:PLUSEND" },
            0.001, 0.01, "ALL", 100.0
        });
    };

    // The interactive session
    //-------------------------------------------------------------------------

    LOG(INFO) << "You are now in the interactive configuration session";
    LOG(INFO) << "Note that this only provides the very basic configuration of the system."
        " For more detailed configuration, read the documents and use the configuration file.";

    LOG(NOTE) << "##################################################";
    LOG(NOTE) << " Units in MEDYAN are nm, second, pN, and pN*nm";
    LOG(NOTE) << "##################################################";

    //---------- basics ----------
    cout << endl;
    LOG(STEP) << "Basic simulation parameters";
    LOG(INFO) << "The geometry information";
    LOG(INFO) << "Note: Network size = compartment size (500nm by default) .* (NX, NY, NZ)";
    cout << "Num compartments in x direction (default 2): "; conf.geoParams.NX = getInt(2);
    cout << "Num compartments in y direction (default 2): "; conf.geoParams.NY = getInt(2);
    cout << "Num compartments in z direction (default 2): "; conf.geoParams.NZ = getInt(2);

    cout << "Boundary shape {CUBIC (default), SPHERICAL, CYLINDER}: ";
    conf.boundParams.boundaryType.boundaryShape = chooseStringAmong({ "CUBIC", "SPHERICAL", "CYLINDER" });
    if(
        const auto& s = conf.boundParams.boundaryType.boundaryShape;
        s == "SPHERICAL" || s == "CYLINDER"
    ) {
        cout << "Diameter of boundary (default 1000): "; conf.boundParams.diameter = getDouble(1000.0);
    }

    cout << endl;
    LOG(INFO) << "Run time specification";
    cout << "Total running time (default 1000): "; conf.chemParams.chemistryAlgorithm.runTime = getDouble(1000.0);
    cout << "Snapshot output time interval (default 1.0): "; conf.chemParams.chemistryAlgorithm.snapshotTime = getDouble(1.0);

    cout << endl;
    LOG(INFO) << "Initial filaments (applicable to type 0 filament only)";
    cout << "Number of filaments (default 30): "; conf.filamentSetup.numFilaments = getInt(30);
    cout << "Number of cylinders of each filament (default 1): "; conf.filamentSetup.filamentLength = getInt(1);

    //---------- chemistry ----------
    cout << endl;
    LOG(STEP) << "Basic chemistry parameters";
    cout << "Initial diffusing actin copy number (default 5000): ";
    conf.chemistryData.speciesDiffusing.push_back({
        "A", getInt(5000), 80.0, 0.0, 0.0, "REG",
        0, "NONE", 0.0
    });
    cout << "Enable filament polymerization/depolymerization (y/n) (default y): ";
    if(chooseYesNo(true)) {
        conf.chemistryData.polymerizationReactions.push_back({
            { "A:DIFFUSING", "PA:PLUSEND" },
            { "AF:FILAMENT", "PA:PLUSEND" },
            0.154, 0.0, "NA"
        });
        conf.chemistryData.polymerizationReactions.push_back({
            { "A:DIFFUSING", "MA:MINUSEND" },
            { "AF:FILAMENT", "MA:MINUSEND" },
            0.017, 0.0, "NA"
        });
        conf.chemistryData.depolymerizationReactions.push_back({
            { "AF:FILAMENT", "PA:PLUSEND" },
            { "A:DIFFUSING", "PA:PLUSEND" },
            1.4, 0.0, "NA"
        });
        conf.chemistryData.depolymerizationReactions.push_back({
            { "AF:FILAMENT", "MA:MINUSEND" },
            { "A:DIFFUSING", "MA:MINUSEND" },
            0.8, 0.0, "NA"
        });
    }

    cout << endl;
    LOG(INFO) << "  Myosin - Non-muscle mysoin IIA";
    LOG(INFO) << "  Linker - alpha-actinin crosslinker";
    LOG(INFO) << "  Bracher - Arp2/3 brancher";
    cout << "Add myosin (y/n) (default y): "; shouldAddMyosin = chooseYesNo(true);
    if(shouldAddMyosin) {
        cout << "Initial diffusing myosin copy number (default 50): ";
        conf.chemistryData.speciesDiffusing.push_back({
            "MD", getInt(50), 0.8, 0.0, 0.0, "REG",
            0, "NONE", 0.0
        });
    }
    cout << "Add linker (y/n) (default y): "; shouldAddLinker = chooseYesNo(true);
    if(shouldAddLinker) {
        cout << "Initial diffusing linker copy number (default 500): ";
        conf.chemistryData.speciesDiffusing.push_back({
            "LD", getInt(500), 8.0, 0.0, 0.0, "REG",
            0, "NONE", 0.0
        });
    }
    cout << "Add brancher (y/n) (default y): "; shouldAddBrancher = chooseYesNo(true);
    if(shouldAddBrancher) {
        cout << "Initial diffusing brancher copy number (default 20): ";
        conf.chemistryData.speciesDiffusing.push_back({
            "BD", getInt(20), 80.0, 1.0, 0.0, "REG",
            0, "NONE", 0.0
        });
    }

    //---------- force field ----------
    cout << endl;
    LOG(STEP) << "Force field parameters";
    LOG(INFO) << "Energy minimization time interval:";
    LOG(INFO) << "  Use a lower value if: 1. simulation fails or generates warnings";
    LOG(INFO) << "                        2. has very fast chemical reactions";
    LOG(INFO) << "  Recommend value: 0.001 - 0.05";
    cout << "Energy minimization interval (default 0.01): ";
    conf.chemParams.chemistryAlgorithm.neighborListTime =
        conf.chemParams.chemistryAlgorithm.minimizationTime = getDouble(0.01);

    cout << endl;
    LOG(INFO) << "Actin force fields (stretching and bending are already enabled)";
    cout << "Enable volume exclusion (y/n) (default y): ";
    if(chooseYesNo(true)) {
        conf.mechParams.mechanicsFFType.VolumeFFType = "REPULSION";
        conf.mechParams.VolumeCutoff = 108.0;
        conf.mechParams.VolumeK = {1e5};
    }

    //---------- dynamic rates ----------
    cout << endl;
    LOG(STEP) << "Dynamic rate parameters";
    LOG(INFO) << "Dynamic rate model - Brownian Ratchet";
    cout << "Enable Brownian Ratchet (y/n) (default y): ";
    if(chooseYesNo(true)) {
        conf.dyRateParams.dynamicRateType.dFPolymerizationType = "BROWRATCHET";
        conf.dyRateParams.dFilPolymerizationCharLength = {2.7};
    }

    if(shouldAddMyosin) {
        cout << "Enable motor catch bond (y/n) (default y): ";
        if(chooseYesNo(true)) {
            conf.dyRateParams.dynamicRateType.dMUnbindingType = "LOWDUTYCATCH";
            conf.dyRateParams.dMotorUnbindingCharForce = {12.62};
        }
        cout << "Enable motor walking stall (y/n) (default y): ";
        if(chooseYesNo(true)) {
            conf.dyRateParams.dynamicRateType.dMWalkingType = "LOWDUTYSTALL";
            conf.dyRateParams.dMotorWalkingCharForce = {90.0};
        }
    }
    if(shouldAddLinker) {
        cout << "Enable linker slip bond (y/n) (default y): ";
        if(chooseYesNo(true)) {
            conf.dyRateParams.dynamicRateType.dLUnbindingType = "SLIP";
            conf.dyRateParams.dLinkerUnbindingCharLength = {0.24};
        }
    }

    //---------- post processing ----------
    if(shouldAddMyosin) addMyosin();
    if(shouldAddLinker) addLinker();
    if(shouldAddBrancher) addBrancher();

    cout << endl;
    LOG(STEP) << "Output files settings";
    cout << "The directory for generated files: "; fileDirectory = getString("");
    if(fileDirectory.empty()) {
        fileDirectory = filesystem::current_path();
        LOG(INFO) << "Using current directory: " << fileDirectory;
    }
    conf.metaParams.inputDirectory = fileDirectory;
    cout << "The name of the system input file: "; conf.metaParams.systemInputFile = fileDirectory / getString("");
    cout << "The name of the chemistry input file: "; conf.chemParams.chemistrySetup.inputFile = getString("");

    // File generation
    //-------------------------------------------------------------------------
    cout << endl;
    LOG(STEP) << "Configuration complete! Generating files...";
    helper.generateInput(conf);

}

} // namespace medyan

#endif
