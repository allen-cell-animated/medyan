#ifndef MEDYAN_MedyanArgs_hpp
#define MEDYAN_MedyanArgs_hpp

#include <string>
#include <thread>

#include "Rand.h"
#include "Util/Io/CmdParse.hpp"
#include "Util/Io/Log.hpp"
#include "utility.h" // rdtsc

struct MedyanCmdInitResult {

    // File and directories
    std::string inputFile;
    std::string inputDirectory;
    std::string outputDirectory;

    // Threadings
    int numThreads = -1;

    // RNG
    bool rngSeedFixed = false;
    unsigned long long rngSeed = 0;

};

inline auto medyanInitFromCommandLine(int argc, char** argv) {
    using namespace cmdparse;

    // The initialization result output
    MedyanCmdInitResult res;

    // Parsing the command line
    {
        Command cmdMain("MEDYAN", "");

        cmdMain.addOptionWithVar('s', "", "file", "System input file", true, res.inputFile);
        cmdMain.addOptionWithVar('i', "", "path", "Input directory", true, res.inputDirectory);
        cmdMain.addOptionWithVar('o', "", "path", "Output directory", true, res.outputDirectory);
        cmdMain.addOption(0, "seed-fixed", "seed", "Fixed random generator seed", false,
            [&](const std::string& arg) {
                res.rngSeedFixed = true;
                VariableWrite<unsigned long long>{std::string("seed")}(res.rngSeed, arg);
            }
        );
        cmdMain.addOptionWithVar('t', "", "int", "Thread count (0 for auto)", false, res.numThreads);
        cmdMain.addHelp();

        try {
            cmdMain.parse(argc, argv);
        } catch (const CommandLogicError& e) {
            std::cerr << e.what() << std::endl;
            // Internal error, no help message generated.
            throw;
        } catch (const ParsingError& e) {
            std::cerr << e.what() << std::endl;
            cmdMain.printUsage();
            throw;
        } catch (const ValidationError& e) {
            std::cerr << e.what() << std::endl;
            cmdMain.printUsage();
            throw;
        }
    }

    //-------------------------------------------------------------------------
    // Post processing
    //-------------------------------------------------------------------------

    // Initialize the logger
    //---------------------------------
    ::medyan::logger::Logger::defaultLoggerInitialization();
    // TODO remove logger initialization here

    // Number of threads
    //---------------------------------
    // The actual number of threads will also be affected by the result of
    // std::thread::hardware_concurrency()
    //
    // Let N1 be the user input, N2 be the result of calling hardward_concurrency,
    // and N be the final number of threads. Then
    //
    //   - N1 < 0: (Treat as unspecified) N = 1
    //   - N1 == 0
    //       - N2 == 0: N = 1; Show warning
    //       - N2 > 0: N = N2
    //   - N1 > 0
    //       - N2 == 0: N = N1
    //       - 0 < N2 < N1: N = N1; Show warning
    //       - N2 >= N1: N = N1
    if(res.numThreads < 0) {
        res.numThreads = 1;
    } else {
        unsigned hc = std::thread::hardware_concurrency();

        if(res.numThreads == 0) {
            if(hc == 0) {
                res.numThreads = 1;
                LOG(WARNING) << "Cannot auto determine the number of threads.";
            } else {
                res.numThreads = hc;
            }
        } else {
            if(hc > 0 && hc < res.numThreads) {
                LOG(WARNING) << "The number of threads specified is greater than hardware concurrency (" << hc << ").";
            }
        }
    }
    LOG(INFO) << "The program will utilize " << res.numThreads << " thread(s) for computing.";

    // Seed global random generator
    //---------------------------------
    if(!res.rngSeedFixed) {
        res.rngSeed = rdtsc();
    }
    LOG(DEBUG) << "Global RNG seed: " << res.rngSeed;
    Rand::eng.seed(res.rngSeed);


    return res;

}

#endif
