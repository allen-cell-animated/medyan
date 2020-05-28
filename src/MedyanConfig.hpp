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
*/

#include <filesystem>
#include <stdexcept>

#include "Parser.h"
#include "SysParams.h"
#include "Util/Io/Log.hpp"

namespace medyan {

// Definition of simulation configuration
struct SimulConfig {

    // The MetaParams class, used to store the source of config file, etc
    struct MetaParams {
        std::filesystem::path systemInputFile;

        // Directory of other input files
        std::filesystem::path inputDirectory;
    };

    MetaParams     metaParams;

    // Parameters from the system input
    MechParams     mechParams;
    ChemParams     chemParams;
    GeoParams      geoParams;
    BoundParams    boundParams;
    DyRateParams   dyRateParams;
    SpecialParams  specialParams;
    BubbleSetup    bubbleSetup;
    FilamentSetup  filamentSetup;

    // Parameters from other inputs
    ChemistryData  chemistryData;
    BubbleData     bubbleData;
    FilamentData   filamentData;

    // Read system configuration from input
    void readSystemConfig(std::istream& is) {
        geoParams      = SystemParser::readGeoParams(is);
        boundParams    = SystemParser::readBoundParams(is, geoParams);
        mechParams     = SystemParser::readMechParams(is);
        chemParams     = SystemParser::readChemParams(is, geoParams);
        dyRateParams   = SystemParser::readDyRateParams(is);
        bubbleSetup    = SystemParser::readBubbleSetup(is);
        filamentSetup  = SystemParser::readFilamentSetup(is);
        specialParams  = SystemParser::readSpecialParams(is);
    }

    // Read chemistry input
    void readChemistryConfig(std::istream& is) {
        chemistryData = ChemistryParser::readChemistryInput(is, chemParams);
    }

    // Read bubble input
    void readBubbleConfig(std::istream& is) {
        bubbleData = BubbleParser::readBubbles(is);
    }

    // Read filament input
    void readFilamentConfig(std::istream& is) {
        filamentData = FilamentParser::readFilaments(is);
    }
};

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
        ReadFile file(systemInputFile);
        conf.readSystemConfig(file.ifs);
    }

    // Read chemistry input
    if(conf.chemParams.chemistrySetup.inputFile.empty()) {
        LOG(FATAL) << "Need to specify a chemical input file. Exiting.";
        throw std::runtime_error("No chemistry input file specified.");
    }
    else {
        ReadFile file(inputDirectory / conf.chemParams.chemistrySetup.inputFile);
        conf.readChemistryConfig(file.ifs);
    }

    // Read bubble input
    if(!conf.bubbleSetup.inputFile.empty()) {
        ReadFile file(inputDirectory / conf.bubbleSetup.inputFile);
        conf.readBubbleConfig(file.ifs);
    }

    // Read filament input
    if(!conf.filamentSetup.inputFile.empty()) {
        ReadFile file(inputDirectory / conf.filamentSetup.inputFile);
        conf.readFilamentConfig(file.ifs);
    }

    return conf;
}

} // namespace medyan

#endif
