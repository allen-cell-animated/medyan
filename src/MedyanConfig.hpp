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


  Functions
  ---------
    SimulConfig readSimulConfig(systemInputFile, inputDirectory)

      Read from configuration file recursively.

*/

#include <filesystem>
#include <stdexcept>

#include "Parser.h"
#include "SysParams.h"
#include "Util/Io/Log.hpp"

namespace medyan {


//---------------------------------
// Read from input
//---------------------------------

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


// Read chemistry input
inline void readChemistryConfig(SimulConfig& sc, std::istream& is) {
    sc.chemistryData = ChemistryParser::readChemistryInput(is, sc.chemParams);
}

// Read bubble input
inline void readBubbleConfig(SimulConfig& sc, std::istream& is) {
    sc.bubbleData = BubbleParser::readBubbles(is);
}

// Read filament input
inline void readFilamentConfig(SimulConfig& sc, std::istream& is) {
    sc.filamentData = FilamentParser::readFilaments(is);
}

// Read simulation configuration from file
inline SimulConfig readSimulConfig(
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
        SystemParser sp;
        sp.parseSystemInput(conf, readFileToString(systemInputFile));
    }

    // Read chemistry input
    if(conf.chemParams.chemistrySetup.inputFile.empty()) {
        LOG(FATAL) << "Need to specify a chemical input file. Exiting.";
        throw std::runtime_error("No chemistry input file specified.");
    }
    else {
        ReadFile file(inputDirectory / conf.chemParams.chemistrySetup.inputFile);
        readChemistryConfig(conf, file.ifs);
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

} // namespace medyan

#endif
