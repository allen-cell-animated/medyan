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

} // namespace medyan

#endif
