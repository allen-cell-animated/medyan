#include "core/io/command_line.h"

#include "core/globals.h"
#include "core/io/command_line_opt.h"
#include "core/io/log.h"

namespace medyan {
namespace commandline {

/// MEDYAN read from command line input
void initializeFromCommandLine(int argc, char** argv) {

    /**************************************************************************
    Parsing and reading from command line input
    **************************************************************************/

    // Preset global variables
    Global::global().mode = GlobalVar::RunMode::Simulation;

    bool runHelp = false;

    // Normal run options
    Option1<std::string> opInputFile {"System input file name", "-s", "SYS_INPUT", &Global::global().systemInputFileName};
    Option1<std::string> opInputDir {"Input directory", "-i", "INPUT_DIR", &Global::global().inputDirectory};
    Option1<std::string> opOutputDir {"Output directory", "-o", "OUTPUT_DIR", &Global::global().outputDirectory};
    opInputFile.require();
    opInputDir.require();
    opOutputDir.require();
    Option1<std::string> opLogFile {"Name of log file", "-l,--log", "LOG_FILE", &Global::global().logFileName};

    Command cmdAnalyze {"Run analysis instead of simulation", "analyze", nullptr,
        []()->bool { Global::global().mode = GlobalVar::RunMode::Analysis; return true; }};

    PosHolder holderRun({&opInputFile, &opInputDir, &opOutputDir, &opLogFile}, {&cmdAnalyze}); holderRun.require();

    // Help options
    Option0 opHelp{ "Print help message", "-h,--help", &runHelp }; opHelp.require();
    PosHolder holderHelp({&opHelp}, {}); holderHelp.require();

    // Main command
    PosMutuallyExclusive meMain({ &holderRun, &holderHelp }); meMain.require();
    Command cmdMain {"", "MEDYAN", &meMain, []{return true;}}; cmdMain.setMain();

    // Main parsing
    if(!commandLineParse(argc, argv, cmdMain)) {
        cmdMain.printUsage();
        exit(EXIT_FAILURE);
    }

    /**************************************************************************
    Initializations
    **************************************************************************/

    if(runHelp) {
        cmdMain.printUsage();
        exit(EXIT_SUCCESS);
    }

    // Initialize logger
    MEDYAN_LOG_DEFAULT_CONFIGURATION(Global::readGlobal().outputDirectory + "/" + Global::readGlobal().logFileName);

}


} // namespace commandline
} // namespace medyan
