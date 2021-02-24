#ifndef MEDYAN_MedyanArgs_hpp
#define MEDYAN_MedyanArgs_hpp

#include <string>
#include <thread>

#include "Core/Globals.hpp"
#include "Rand.h"
#include "Util/Io/CmdParse.hpp"
#include "Util/Io/Log.hpp"
#include "utility.h" // rdtsc

enum class MedyanRunMode {
    simulation,
    gui,
    analyze,
    config,
    test
};

struct MedyanCmdInitResult {

    // Remaining arguments
    int argpNext = 0; // If the parsing is complete, this points to invalid location

    // Run mode
    MedyanRunMode runMode = MedyanRunMode::simulation;

    // File and directories
    std::string inputFile;
    std::string inputDirectory;
    std::string outputDirectory;

    // GUI
    bool guiEnabled = false;

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

        cmdMain.addOption(makeOptionWithVar('s', "", "file", "System input file", false, res.inputFile));
        cmdMain.addOption(makeOptionWithVar('i', "", "path", "Input directory", false, res.inputDirectory));
        cmdMain.addOption(makeOptionWithVar('o', "", "path", "Output directory", false, res.outputDirectory));
        cmdMain.addOption(Option(0, "seed-fixed", "seed", "Fixed random generator seed", false,
            [&](const Command&, const std::string& arg) {
                res.rngSeedFixed = true;
                VariableWrite<unsigned long long>{std::string("seed")}(res.rngSeed, arg);
            }
        ));
        cmdMain.addOption(Option(0, "gui", "Enable GUI", false,
            [&](const Command&) {
                res.guiEnabled = true;
            }
        ));
        cmdMain.addOption(makeOptionWithVar('t', "", "int", "Thread count (0 for auto)", false, res.numThreads));
        cmdMain.addHelp();

        // Add gui command
        Command& cmdGui = cmdMain.addCommand(
            "gui", "GUI mode",
            [&] { res.runMode = MedyanRunMode::gui; }
        );
        cmdGui.addHelp();

        // Add analyze command
        Command& cmdAnalyze = cmdMain.addCommand(
            "analyze", "Analyze simulation output",
            [&] { res.runMode = MedyanRunMode::analyze; }
        );
        cmdAnalyze.addOption(Option(0, "bond-frame", "frame", "Frame of membrane topology information", false,
            [&](const Command&, const std::string& arg) {
                if(arg == "all") {
                    medyan::globalMutable().analyzeMembraneBondAllFrames = true;
                } else {
                    VariableWrite< size_t >{"bond-frame"}(medyan::globalMutable().analyzeMembraneBondFrame, arg);
                }
            }
        ));
        cmdAnalyze.addOption(makeOptionWithVar(0, "frame-interval", "int", "Interval of frames", false, medyan::globalMutable().analyzeFrameInterval));
        cmdAnalyze.addHelp();

        // Add interactive configuration command
        auto& cmdConfig = cmdMain.addCommand(
            "config", "Interactive system configuration.",
            [&] { res.runMode = MedyanRunMode::config; }
        );

        // Add test command
        auto& cmdTest = cmdMain.addCommand(
            "test", "Run MEDYAN tests.",
            [&] { res.runMode = MedyanRunMode::test; }
        );
        cmdTest.setTerminating(true);

        // Add validation
        cmdMain.setValidation([&] {
            if(res.runMode == MedyanRunMode::simulation) {
                if(res.inputFile.empty()) {
                    throw ValidationError("Must specify the system input file (-s)");
                }
                if(res.inputDirectory.empty()) {
                    throw ValidationError("Must specify the input directory (-i)");
                }
                if(res.outputDirectory.empty()) {
                    throw ValidationError("Must specify the output directory (-o)");
                }
            }
        });

        try {

            cmdMain.ruleCheck();
            res.argpNext = cmdMain.parse(argc, argv);

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
