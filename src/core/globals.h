#ifndef MEDYAN_CORE_GLOBALS_H
#define MEDYAN_CORE_GLOBALS_H

#include <string>

namespace medyan {
    struct GlobalVar {
        enum class RunMode {
            Simulation,
            Analysis
        } mode; ///< The core mode at which the program runs

        std::string systemInputFile; // Path + name
        std::string inputDirectory;
        std::string outputDirectory;

        std::string logFileName {"medyan.log"};

        bool randomGenSeedFixed = false;
        long long randomGenSeed; // Only used in fixed random seed
    
    };

    /// GlobalVar manager
    class Global {
        static GlobalVar _global;
    public:
        static const GlobalVar& readGlobal() { return _global; }
        static GlobalVar& global() { return _global; }
    };
}


#endif