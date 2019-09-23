#ifndef MEDYAN_Core_Globals_Hpp
#define MEDYAN_Core_Globals_Hpp

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
    unsigned long long randomGenSeed; // The global random seed

    // Analyze specific
    std::size_t analyzeMembraneBondFrame = 0;
    int analyzeFrameInterval = 1; // must be >= 1

};

inline auto& globalMutable() { static GlobalVar g; return g; }
inline const auto& global() { return globalMutable(); }

} // namespace medyan


#endif