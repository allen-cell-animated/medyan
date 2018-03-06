#ifndef MEDYAN_CORE_GLOBALS_H
#define MEDYAN_CORE_GLOBALS_H

#include <string>

namespace medyan {
    constexpr const char* programName = "MEDYAN";
    
    struct GlobalVar {
        /**
         * The core mode at which the program runs
         * 
         * 0: Simulation
         * 1: Analyzing
         */
        int mode;

        std::string systemInputFileName;
        std::string inputDirectory;
        std::string outputDirectory;
    
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